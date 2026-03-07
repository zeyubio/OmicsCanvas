#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OmicsCanvas - Genome-wide Circos-like plot (bins-based, matplotlib only)
=======================================================================

Reset requirements
------------------
1) Y-axis ticks: only "0" and "max", drawn ONLY once on the first chromosome.
2) Genes: draw uniform green ticks at gene midpoints on the innermost ideogram line.
3) Legend: upper-right.
4) Chromosome names: inside the ideogram.

Dependencies
------------
python >=3.8
pip install numpy pandas matplotlib pysam
"""

from __future__ import annotations

import os
import re
import math
import argparse
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42   # TrueType
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['svg.fonttype'] = 'none'  # SVG 中保留文字
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 14


# -----------------------------
# Utils
# -----------------------------
def parse_list(arg: Optional[str]) -> List[str]:
    if not arg:
        return []
    return [x.strip() for x in arg.split(",") if x.strip()]


def resolve_paths_with_dir(paths: List[str], base_dir: Optional[str]) -> List[str]:
    if not paths:
        return []
    bd = os.path.expanduser(base_dir) if base_dir else None
    out: List[str] = []
    for p in paths:
        p2 = os.path.expanduser(p)
        if os.path.isabs(p2) or p2.startswith(("s3://", "gs://")):
            out.append(p2)
        else:
            out.append(os.path.join(bd, p2) if bd else p2)
    return out


def fmt_axis_value(v: float) -> str:
    try:
        v = float(v)
    except Exception:
        return str(v)
    if not np.isfinite(v):
        return "NA"
    av = abs(v)
    if av >= 1e9:
        return f"{v/1e9:.2g}G"
    if av >= 1e6:
        return f"{v/1e6:.2g}M"
    if av >= 1e3:
        return f"{v/1e3:.2g}k"
    if av >= 1:
        return f"{v:.3g}"
    if av >= 1e-3:
        return f"{v:.3g}"
    return f"{v:.2e}"


def _default_color_cycle() -> List[str]:
    return [
        "#ee6666", "#5470C6", "#91CC75", "#FAC858", "#73C0DE",
        "#3BA272", "#FC8452", "#9A60B4", "#EA7CCC", "#2E91E5",
    ]


# -----------------------------
# Chrom sizes
# -----------------------------
def read_chrom_sizes_from_whitespace_2col(path: str) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = re.split(r"\s+", s)
            if len(parts) < 2:
                continue
            try:
                sizes[parts[0]] = int(parts[1])
            except ValueError:
                continue
    return sizes


def infer_chrom_sizes_from_bam(bam_path: str) -> Dict[str, int]:
    import pysam
    sizes: Dict[str, int] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, ln in zip(bam.references, bam.lengths):
            sizes[str(chrom)] = int(ln)
    return sizes


# -----------------------------
# Chrom selection + layout
# -----------------------------
def _starts_with_prefix(name: str, prefix: str, ignore_case: bool) -> bool:
    if not prefix:
        return True
    return name.lower().startswith(prefix.lower()) if ignore_case else name.startswith(prefix)


def pick_chromosomes(
    chrom_sizes: Dict[str, int],
    chroms: Optional[List[str]] = None,
    chrom_prefix: str = "Chr",
    prefix_ignore_case: bool = False,
    exclude_keywords: str = "random,Un,un,scaffold,Scaffold,chrUn,chrM,MT,chloroplast,plastid,mitochond",
    max_chroms: int = 25,
    min_len: int = 1,
) -> List[str]:
    if chroms:
        chosen = [c for c in chroms if c in chrom_sizes]
        if not chosen:
            raise ValueError("None of the requested --chroms exist in chrom sizes.")
        return chosen

    kws = [k.strip() for k in exclude_keywords.split(",") if k.strip()]
    ere = re.compile("|".join([re.escape(k) for k in kws]), flags=re.IGNORECASE) if kws else None

    items: List[Tuple[str, int]] = []
    for c, ln in chrom_sizes.items():
        if ln < min_len:
            continue
        if ere and ere.search(c):
            continue
        if not _starts_with_prefix(c, chrom_prefix, prefix_ignore_case):
            continue
        items.append((c, ln))

    if not items:
        tmp: List[Tuple[str, int]] = []
        for c, ln in chrom_sizes.items():
            if ln < min_len:
                continue
            if ere and ere.search(c):
                continue
            tmp.append((c, ln))
        items = tmp if tmp else list(chrom_sizes.items())

    items = sorted(items, key=lambda x: x[1], reverse=True)
    return [c for c, _ in items[:max_chroms]]


def build_layout_units(
    chrom_list: List[str],
    chrom_sizes: Dict[str, int],
    bin_size: int,
    gap_frac: float = 0.03,
    gap_bins: Optional[int] = None,
    gap_after_last: bool = True,
) -> Tuple[Dict[str, int], Dict[str, int], int]:
    chrom_start: Dict[str, int] = {}
    chrom_nbins: Dict[str, int] = {}
    cur = 0

    def _gap_for(nb: int) -> int:
        return int(max(1, round(nb * gap_frac))) if gap_bins is None else int(gap_bins)

    for i, c in enumerate(chrom_list):
        nb = int(math.ceil(chrom_sizes[c] / bin_size))
        chrom_start[c] = cur
        chrom_nbins[c] = nb
        cur += nb
        if i != len(chrom_list) - 1:
            cur += _gap_for(nb)

    if gap_after_last and chrom_list:
        nb_last = chrom_nbins[chrom_list[-1]]
        cur += _gap_for(nb_last)

    return chrom_start, chrom_nbins, cur


# -----------------------------
# Data
# -----------------------------
def read_gene_midpoints(gene_bed: str, chrom_list: List[str]) -> Dict[str, np.ndarray]:
    df = pd.read_csv(
        gene_bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": "string", "start": "int64", "end": "int64"},
        comment="#",
    )
    df = df[df["chrom"].isin(chrom_list)].copy()
    out: Dict[str, np.ndarray] = {c: np.array([], dtype=np.int64) for c in chrom_list}
    if df.empty:
        return out
    mid = ((df["start"].values + df["end"].values) // 2).astype(np.int64)
    chrom = df["chrom"].astype(str).values
    for c in chrom_list:
        out[c] = mid[chrom == c]
    return out


def bin_methylation_cx(
    cx_path: str,
    chrom_list: List[str],
    chrom_sizes: Dict[str, int],
    bin_size: int,
    min_bin_depth: int = 1,
    chunksize: int = 2_000_000,
) -> Dict[str, np.ndarray]:
    meth_sum = {c: np.zeros(int(math.ceil(chrom_sizes[c] / bin_size)), dtype=np.float64) for c in chrom_list}
    dep_sum  = {c: np.zeros(int(math.ceil(chrom_sizes[c] / bin_size)), dtype=np.float64) for c in chrom_list}

    colnames = ["chrom", "pos", "meth", "depth"]
    reader = pd.read_csv(
        cx_path,
        sep="\t",
        header=None,
        names=colnames,
        usecols=[0, 1, 2, 3],
        dtype={0: "string", 1: "int64", 2: "int64", 3: "int64"},
        comment="#",
        chunksize=chunksize,
        engine="c",
    )
    chrom_set = set(chrom_list)

    for chunk in reader:
        chunk = chunk.dropna()
        chunk = chunk[chunk["depth"] > 0]
        if chunk.empty:
            continue
        chunk = chunk[chunk["chrom"].isin(chrom_set)]
        if chunk.empty:
            continue

        bins = ((chunk["pos"].values - 1) // bin_size).astype(np.int64)
        chroms = chunk["chrom"].astype(str).values
        meths = chunk["meth"].values.astype(np.float64)
        deps  = chunk["depth"].values.astype(np.float64)

        for c, b, me, de in zip(chroms, bins, meths, deps):
            if b < meth_sum[c].shape[0]:
                meth_sum[c][b] += me
                dep_sum[c][b]  += de

    ratio: Dict[str, np.ndarray] = {}
    for c in chrom_list:
        d = dep_sum[c]
        me = meth_sum[c]
        r = np.zeros_like(d, dtype=np.float64)
        ok = d >= float(min_bin_depth)
        r[ok] = me[ok] / d[ok]
        ratio[c] = r
    return ratio


def bin_bam_counts(
    bam_path: str,
    chrom_list: List[str],
    chrom_sizes: Dict[str, int],
    bin_size: int,
    pe_count: str = "read1",
    include_duplicates: bool = False,
    mapq: int = 1,
    norm: str = "rpm",
) -> Tuple[Dict[str, np.ndarray], float]:
    import pysam
    out = {c: np.zeros(int(math.ceil(chrom_sizes[c] / bin_size)), dtype=np.float64) for c in chrom_list}

    total = 0.0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for c in chrom_list:
            try:
                it = bam.fetch(c)
            except ValueError:
                continue
            for r in it:
                if r.is_unmapped:
                    continue
                if r.is_secondary or r.is_supplementary:
                    continue
                if r.is_qcfail:
                    continue
                if (not include_duplicates) and r.is_duplicate:
                    continue
                if r.mapping_quality < mapq:
                    continue
                if r.is_paired and pe_count == "read1" and (not r.is_read1):
                    continue
                pos = r.reference_start
                if pos < 0:
                    continue
                b = int(pos // bin_size)
                arr = out.get(c)
                if arr is None or b >= arr.shape[0]:
                    continue
                arr[b] += 1.0
                total += 1.0

    if norm == "rpm" and total > 0:
        scale = 1e6 / total
        for c in chrom_list:
            out[c] *= scale
    return out, total


def scale_values(values: np.ndarray, method: str = "p99", eps: float = 1e-12) -> float:
    v = values[np.isfinite(values)]
    v = v[v >= 0]
    if v.size == 0:
        return eps
    if method == "none":
        return max(float(v.max()), eps)
    if method == "max":
        s = float(v.max())
    elif method == "p95":
        s = float(np.percentile(v, 95))
    else:
        s = float(np.percentile(v, 99))
    if s < eps:
        s = float(v.max())
    return max(s, eps)


# -----------------------------
# Plotting primitives
# -----------------------------
def _plot_track_spokes(ax, chrom_list, chrom_start, chrom_nbins, rotate, one, r0, height, values_by_chrom, clip_val, color, seg_lw):
    from matplotlib.collections import LineCollection
    segs = []
    for c in chrom_list:
        arr = values_by_chrom.get(c)
        if arr is None or arr.size == 0:
            continue
        nb = chrom_nbins[c]
        arr = arr[:nb].astype(np.float64, copy=False)
        u = chrom_start[c] + (np.arange(nb, dtype=np.float64) + 0.5)
        theta = rotate - one * u
        v = np.clip(arr, 0, None)
        denom = clip_val if (clip_val and clip_val > 0) else 1.0
        if clip_val and clip_val > 0:
            v = np.minimum(v, clip_val)
        r2 = r0 + (v / denom) * height
        x1 = np.cos(theta) * r0
        y1 = np.sin(theta) * r0
        x2 = np.cos(theta) * r2
        y2 = np.sin(theta) * r2
        segs.append(np.stack([np.stack([x1, y1], axis=1), np.stack([x2, y2], axis=1)], axis=1))
    if segs:
        segs = np.concatenate(segs, axis=0)
        ax.add_collection(LineCollection(segs, colors=color, linewidths=seg_lw))


def _plot_track_fill(ax, chrom_list, chrom_start, chrom_nbins, rotate, one, r0, height, values_by_chrom, clip_val, color, alpha):
    from matplotlib.collections import PolyCollection
    polys = []
    for c in chrom_list:
        arr = values_by_chrom.get(c)
        if arr is None or arr.size == 0:
            continue
        nb = chrom_nbins[c]
        arr = arr[:nb].astype(np.float64, copy=False)
        v = np.clip(arr, 0, None)
        denom = clip_val if (clip_val and clip_val > 0) else 1.0
        if clip_val and clip_val > 0:
            v = np.minimum(v, clip_val)

        u_edges = chrom_start[c] + np.arange(nb + 1, dtype=np.float64)
        theta_edges = rotate - one * u_edges
        outer_edges = r0 + (np.r_[v, (v[-1] if v.size else 0.0)] / denom) * height
        inner_edges = np.full_like(outer_edges, r0, dtype=np.float64)

        x_out = np.cos(theta_edges) * outer_edges
        y_out = np.sin(theta_edges) * outer_edges
        x_in  = np.cos(theta_edges[::-1]) * inner_edges[::-1]
        y_in  = np.sin(theta_edges[::-1]) * inner_edges[::-1]
        verts = np.vstack([np.column_stack([x_out, y_out]), np.column_stack([x_in, y_in])])
        polys.append(verts)

    if polys:
        ax.add_collection(PolyCollection(polys, facecolors=color, edgecolors="none", alpha=alpha))


def _plot_gene_ticks_on_ideogram(ax, gene_mid_by_chrom, chrom_list, chrom_start,
                                rotate, one, bin_size, ideogram_inner_r,
                                tick_len, color, lw, stride=1):
    from matplotlib.collections import LineCollection
    segs = []
    stride = max(1, int(stride))
    for c in chrom_list:
        mids = gene_mid_by_chrom.get(c)
        if mids is None or len(mids) == 0:
            continue
        mids = mids[::stride]
        u = chrom_start[c] + (mids.astype(np.float64) / float(bin_size))
        theta = rotate - one * u
        x1 = np.cos(theta) * ideogram_inner_r
        y1 = np.sin(theta) * ideogram_inner_r
        x2 = np.cos(theta) * (ideogram_inner_r - tick_len)
        y2 = np.sin(theta) * (ideogram_inner_r - tick_len)
        segs.append(np.stack([np.stack([x1, y1], axis=1), np.stack([x2, y2], axis=1)], axis=1))
    if segs:
        segs = np.concatenate(segs, axis=0)
        ax.add_collection(LineCollection(segs, colors=color, linewidths=lw))


def _draw_track_y_ticks_on_first_chrom(ax, tracks, first_chrom, chrom_start, chrom_nbins,
                                      rotate, one, base_radius, track_height, track_gap,
                                      ytick_angle_deg=0.0,
                                      fontsize=8, tick_lw=0.8):
    if not first_chrom or first_chrom not in chrom_start:
        return
        # Place all y tick labels at a fixed angle (default: right middle, above).
    # This makes the y tick labels visually consistent and avoids drifting with chromosome length.
    th = math.radians(float(ytick_angle_deg))

    off = 0.012
    tx_off = -math.sin(th) * off
    ty_off =  math.cos(th) * off

    for i, tr in enumerate(tracks):
        r0 = base_radius + i * (track_height + track_gap)
        r1 = r0 + track_height
        vmax = float(tr.get("clip_val", tr.get("scale", 0.0)))

        ax.plot([math.cos(th) * r0, math.cos(th) * r1],
                [math.sin(th) * r0, math.sin(th) * r1],
                color="black", lw=tick_lw)

        x0, y0 = math.cos(th) * (r0 - 0.02), math.sin(th) * (r0 - 0.02)
        x1, y1 = math.cos(th) * (r1 + 0.04), math.sin(th) * (r1 + 0.04)

        ax.text(x0 + tx_off, y0 + ty_off, "0", ha="right", va="center", fontsize=fontsize, color="black")
        ax.text(x1 + tx_off, y1 + ty_off, fmt_axis_value(vmax), ha="right", va="center", fontsize=fontsize, color="black")


def plot_circos(
    tracks: List[dict],
    chrom_list: List[str],
    chrom_start: Dict[str, int],
    chrom_nbins: Dict[str, int],
    total_units: int,
    out_fig: str,
    bin_size: int,
    gene_mid_by_chrom: Optional[Dict[str, np.ndarray]] = None,
    gene_tick_color: str = "#00AA00",
    gene_tick_len: float = 0.035,
    gene_tick_lw: float = 0.6,
    gene_tick_stride: int = 1,
    fig_size: float = 10.0,
    dpi: int = 300,
    ideogram_radius: float = 1.0,
    ideogram_width: float = 0.08,
    chr_label_radius: float = 0.86,
    base_radius: float = 1.10,
    track_height: float = 0.20,
    track_gap: float = 0.13,
    rotate_deg: float = 90.0,
    outline_lw: float = 0.8,
    seg_lw: float = 0.25,
    fill_alpha: float = 0.7,
    legend: bool = True,
    legend_fontsize: int = 10,
    ytick_angle_deg: float = 0.0,
):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.collections import PolyCollection

    rotate = math.radians(rotate_deg)
    one = 2 * math.pi / float(total_units)

    fig = plt.figure(figsize=(fig_size, fig_size))
    ax = fig.add_axes([0.05, 0.05, 0.75, 0.90])
    ax.set_aspect("equal")
    ax.axis("off")

    # ideogram band
    r_outer = ideogram_radius
    r_inner = ideogram_radius - ideogram_width

    polys = []
    for c in chrom_list:
        start_u = chrom_start[c]
        nb = chrom_nbins[c]
        u = np.linspace(start_u, start_u + nb, min(max(300, nb // 2), 2500))
        theta = rotate - one * u
        x_out = np.cos(theta) * r_outer
        y_out = np.sin(theta) * r_outer
        x_in  = np.cos(theta[::-1]) * r_inner
        y_in  = np.sin(theta[::-1]) * r_inner
        verts = np.vstack([np.column_stack([x_out, y_out]), np.column_stack([x_in, y_in])])
        polys.append(verts)
    ax.add_collection(PolyCollection(polys, facecolors="#E7F4D8", edgecolors="none", alpha=1.0))

    # ideogram outlines + boundaries + chr labels inside
    for c in chrom_list:
        start_u = chrom_start[c]
        nb = chrom_nbins[c]
        u = np.linspace(start_u, start_u + nb, min(max(300, nb // 2), 2500))
        theta = rotate - one * u
        ax.plot(np.cos(theta) * r_outer, np.sin(theta) * r_outer, color="black", lw=outline_lw)
        ax.plot(np.cos(theta) * r_inner, np.sin(theta) * r_inner, color="black", lw=outline_lw)
        for uu in (start_u, start_u + nb):
            th = rotate - one * uu
            ax.plot([np.cos(th) * r_inner, np.cos(th) * r_outer],
                    [np.sin(th) * r_inner, np.sin(th) * r_outer],
                    color="black", lw=outline_lw)

        mid_u = start_u + nb / 2.0
        thm = rotate - one * mid_u
        label_r = min(float(chr_label_radius), float(r_inner) - 0.10)
        label_r = max(label_r, 0.20)
        tx, ty = np.cos(thm) * label_r, np.sin(thm) * label_r
        rot = (math.degrees(thm) - 90) % 360
        if 90 < rot < 270:
            rot = (rot + 180) % 360
        ax.text(tx, ty, c, ha="center", va="center", fontsize=11, rotation=rot)

    # gene ticks
    if gene_mid_by_chrom is not None:
        _plot_gene_ticks_on_ideogram(
            ax=ax,
            gene_mid_by_chrom=gene_mid_by_chrom,
            chrom_list=chrom_list,
            chrom_start=chrom_start,
            rotate=rotate,
            one=one,
            bin_size=bin_size,
            ideogram_inner_r=r_inner,
            tick_len=gene_tick_len,
            color=gene_tick_color,
            lw=gene_tick_lw,
            stride=gene_tick_stride,
        )

    # tracks
    handles = []
    labels = []
    for i, tr in enumerate(tracks):
        r0 = base_radius + i * (track_height + track_gap)
        height = track_height
        color = tr["color"]
        name = tr["name"]
        values_by_chrom = tr["values_by_chrom"]
        clip_val = float(tr.get("clip_val", tr.get("scale", 0.0)))
        style = tr.get("style", "spokes")

        # outlines per chromosome
        for c in chrom_list:
            start_u = chrom_start[c]
            nb = chrom_nbins[c]
            u = np.linspace(start_u, start_u + nb, min(max(250, nb // 2), 2500))
            theta = rotate - one * u
            ax.plot(np.cos(theta) * r0, np.sin(theta) * r0, color="black", lw=0.55)
            ax.plot(np.cos(theta) * (r0 + height), np.sin(theta) * (r0 + height), color="black", lw=0.55)
            for uu in (start_u, start_u + nb):
                th = rotate - one * uu
                ax.plot([np.cos(th) * r0, np.cos(th) * (r0 + height)],
                        [np.sin(th) * r0, np.sin(th) * (r0 + height)],
                        color="black", lw=0.55)

        if style == "fill":
            _plot_track_fill(ax, chrom_list, chrom_start, chrom_nbins, rotate, one, r0, height, values_by_chrom, clip_val, color, fill_alpha)
            handles.append(Line2D([0], [0], color=color, lw=8))
        else:
            _plot_track_spokes(ax, chrom_list, chrom_start, chrom_nbins, rotate, one, r0, height, values_by_chrom, clip_val, color, seg_lw)
            handles.append(Line2D([0], [0], color=color, lw=2))
        labels.append(name)

    # y ticks on first chromosome
    if chrom_list:
        _draw_track_y_ticks_on_first_chrom(
            ax=ax,
            tracks=tracks,
            first_chrom=chrom_list[0],
            chrom_start=chrom_start,
            chrom_nbins=chrom_nbins,
            rotate=rotate,
            one=one,
            base_radius=base_radius,
            track_height=track_height,
            track_gap=track_gap,
            fontsize=8,
            tick_lw=0.8,
            ytick_angle_deg=ytick_angle_deg,
        )

    # legend upper-right
    if legend and handles:
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.98, 0.98),
                   frameon=False, fontsize=legend_fontsize, handlelength=2.0, labelspacing=0.6)

    rmax = base_radius + len(tracks) * (track_height + track_gap) + 0.35
    lim = max(rmax, ideogram_radius + 0.25)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    fig.savefig(out_fig, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ap.add_argument("--cx-dir", default=None, help="Base directory for CX files.")
    ap.add_argument("--bam-dir", default=None, help="Base directory for omics BAM files.")
    ap.add_argument("--rna-dir", default=None, help="Base directory for RNA BAM files.")

    ap.add_argument("--gene-bed", default=None, help="Gene BED for uniform ticks on ideogram.")
    ap.add_argument("--gene-tick-color", default="#00AA00", help="Gene tick color.")
    ap.add_argument("--gene-tick-len", type=float, default=0.035, help="Gene tick length (radius units).")
    ap.add_argument("--gene-tick-lw", type=float, default=0.6, help="Gene tick linewidth.")
    ap.add_argument("--gene-tick-stride", type=int, default=1, help="Plot every N-th gene tick (1=all).")

    ap.add_argument("--cx-tracks", default=None, help="Comma-separated CX files.")
    ap.add_argument("--cx-names", default=None, help="Comma-separated CX names.")
    ap.add_argument("--cx-colors", default=None, help="Comma-separated CX colors.")

    ap.add_argument("--bam-tracks", default=None, help="Comma-separated omics BAM files.")
    ap.add_argument("--bam-names", default=None, help="Comma-separated omics names.")
    ap.add_argument("--bam-colors", default=None, help="Comma-separated omics colors.")

    ap.add_argument("--rna-tracks", default=None, help="Comma-separated RNA BAM files.")
    ap.add_argument("--rna-names", default=None, help="Comma-separated RNA names.")
    ap.add_argument("--rna-colors", default=None, help="Comma-separated RNA colors.")

    ap.add_argument("--bin-size", type=int, default=100000, help="Bin size in bp.")
    ap.add_argument("--min-bin-depth", type=int, default=10, help="Min depth per bin for CX ratio.")
    ap.add_argument("--cx-chunksize", type=int, default=2_000_000, help="CX chunksize (rows).")

    ap.add_argument("--chrom-sizes", default=None, help="2-col chrom sizes file (whitespace separated).")
    ap.add_argument("--fai", default=None, help="FASTA index (.fai).")

    ap.add_argument("--chroms", default=None, help="Comma-separated chromosomes to plot.")
    ap.add_argument("--chrom-prefix", default="Chr", help="Chromosome prefix to select.")
    ap.add_argument("--chrom-prefix-ignore-case", action="store_true")
    ap.add_argument("--exclude-keywords", default="random,Un,un,scaffold,Scaffold,chrUn,chrM,MT,chloroplast,plastid,mitochond")
    ap.add_argument("--max-chroms", type=int, default=25)
    ap.add_argument("--min-chrom-len", type=int, default=1)

    ap.add_argument("--pe-count", choices=["read1", "all"], default="read1")
    ap.add_argument("--include-duplicates", action="store_true")
    ap.add_argument("--mapq", type=int, default=1)
    ap.add_argument("--norm", choices=["none", "rpm"], default="rpm")

    ap.add_argument("--cx-style", choices=["spokes", "fill"], default="spokes")
    ap.add_argument("--omics-style", choices=["spokes", "fill"], default="fill")
    ap.add_argument("--rna-style", choices=["spokes", "fill"], default="fill")
    ap.add_argument("--fill-alpha", type=float, default=0.7)

    ap.add_argument("--scale-method", choices=["max", "p95", "p99", "none"], default="p99")
    ap.add_argument("--clip-mode", choices=["scale", "none"], default="scale")
    ap.add_argument("--clip-mult", type=float, default=1.0)

    ap.add_argument("--gap-frac", type=float, default=0.06)
    ap.add_argument("--gap-bins", type=int, default=None)
    ap.add_argument("--no-gap-after-last", action="store_true")

    ap.add_argument("--fig-size", type=float, default=10.0)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--fig-format", choices=["pdf", "png"], default="pdf")
    ap.add_argument("--ideogram-radius", type=float, default=1.0)
    ap.add_argument("--ideogram-width", type=float, default=0.08)
    ap.add_argument("--chr-label-radius", type=float, default=0.70)
    ap.add_argument("--base-radius", type=float, default=1.10)
    ap.add_argument("--track-height", type=float, default=0.20)
    ap.add_argument("--track-gap", type=float, default=0.13)
    ap.add_argument("--rotate-deg", type=float, default=0.0)
    ap.add_argument("--ytick-angle-deg", type=float, default=0.0,
                    help="Angle (deg) to place the y tick labels (0 and max) for all tracks. 0=right, 90=top.")
    ap.add_argument("--outline-lw", type=float, default=0.8)
    ap.add_argument("--seg-lw", type=float, default=0.25)

    ap.add_argument("--no-legend", action="store_true")
    ap.add_argument("--legend-fontsize", type=int, default=10)

    ap.add_argument("--out-prefix", required=True)

    args = ap.parse_args()

    cx_paths  = resolve_paths_with_dir(parse_list(args.cx_tracks), args.cx_dir)
    bam_paths = resolve_paths_with_dir(parse_list(args.bam_tracks), args.bam_dir)
    rna_paths = resolve_paths_with_dir(parse_list(args.rna_tracks), args.rna_dir)

    # chrom sizes
    if args.chrom_sizes:
        chrom_sizes = read_chrom_sizes_from_whitespace_2col(args.chrom_sizes)
    elif args.fai:
        chrom_sizes = read_chrom_sizes_from_whitespace_2col(args.fai)
    else:
        any_bam = (bam_paths + rna_paths)
        if any_bam:
            chrom_sizes = infer_chrom_sizes_from_bam(any_bam[0])
        else:
            raise SystemExit("ERROR: Need --chrom-sizes/--fai or a BAM to infer sizes.")
    if not chrom_sizes:
        raise SystemExit("ERROR: empty chrom sizes (check --fai/--chrom-sizes).")

    chrom_list = pick_chromosomes(
        chrom_sizes,
        chroms=parse_list(args.chroms) if args.chroms else None,
        chrom_prefix=args.chrom_prefix,
        prefix_ignore_case=args.chrom_prefix_ignore_case,
        exclude_keywords=args.exclude_keywords,
        max_chroms=args.max_chroms,
        min_len=args.min_chrom_len,
    )

    chrom_start, chrom_nbins, total_units = build_layout_units(
        chrom_list,
        chrom_sizes,
        args.bin_size,
        gap_frac=args.gap_frac,
        gap_bins=args.gap_bins,
        gap_after_last=(not args.no_gap_after_last),
    )

    def make_clip_val(scale: float) -> float:
        if args.clip_mode == "none":
            return scale
        return max(scale * float(args.clip_mult), 1e-12)

    colors = _default_color_cycle()
    color_i = 0
    def next_color():
        nonlocal color_i
        c = colors[color_i % len(colors)]
        color_i += 1
        return c

    tracks: List[dict] = []

    # CX tracks
    cx_names = parse_list(args.cx_names)
    cx_colors = parse_list(args.cx_colors)
    if cx_paths:
        if cx_names and len(cx_names) != len(cx_paths):
            raise SystemExit("ERROR: --cx-names count must match --cx-tracks.")
        if cx_colors and len(cx_colors) != len(cx_paths):
            raise SystemExit("ERROR: --cx-colors count must match --cx-tracks.")
        for i, fp in enumerate(cx_paths):
            if not os.path.exists(fp):
                raise SystemExit(f"ERROR: CX not found: {fp}")
            name = cx_names[i] if cx_names else os.path.basename(fp)
            color = cx_colors[i] if cx_colors else next_color()
            vals = bin_methylation_cx(fp, chrom_list, chrom_sizes, args.bin_size, min_bin_depth=args.min_bin_depth, chunksize=args.cx_chunksize)
            all_vals = np.concatenate([vals[c] for c in chrom_list]) if chrom_list else np.array([0.0])
            sc = scale_values(all_vals, method=args.scale_method)
            tracks.append({"name": name, "color": color, "values_by_chrom": vals, "scale": sc, "clip_val": make_clip_val(sc), "style": args.cx_style})

    # BAM omics tracks
    bam_names = parse_list(args.bam_names)
    bam_colors = parse_list(args.bam_colors)
    if bam_paths:
        if bam_names and len(bam_names) != len(bam_paths):
            raise SystemExit("ERROR: --bam-names count must match --bam-tracks.")
        if bam_colors and len(bam_colors) != len(bam_paths):
            raise SystemExit("ERROR: --bam-colors count must match --bam-tracks.")
        for i, fp in enumerate(bam_paths):
            if not os.path.exists(fp):
                raise SystemExit(f"ERROR: BAM not found: {fp}")
            name = bam_names[i] if bam_names else os.path.basename(fp)
            color = bam_colors[i] if bam_colors else next_color()
            vals, _ = bin_bam_counts(fp, chrom_list, chrom_sizes, args.bin_size, pe_count=args.pe_count, include_duplicates=args.include_duplicates, mapq=args.mapq, norm=args.norm)
            all_vals = np.concatenate([vals[c] for c in chrom_list]) if chrom_list else np.array([0.0])
            sc = scale_values(all_vals, method=args.scale_method)
            tracks.append({"name": name, "color": color, "values_by_chrom": vals, "scale": sc, "clip_val": make_clip_val(sc), "style": args.omics_style})

    # RNA tracks
    rna_names = parse_list(args.rna_names)
    rna_colors = parse_list(args.rna_colors)
    if rna_paths:
        if rna_names and len(rna_names) != len(rna_paths):
            raise SystemExit("ERROR: --rna-names count must match --rna-tracks.")
        if rna_colors and len(rna_colors) != len(rna_paths):
            raise SystemExit("ERROR: --rna-colors count must match --rna-tracks.")
        for i, fp in enumerate(rna_paths):
            if not os.path.exists(fp):
                raise SystemExit(f"ERROR: RNA BAM not found: {fp}")
            name = rna_names[i] if rna_names else os.path.basename(fp)
            color = rna_colors[i] if rna_colors else next_color()
            vals, _ = bin_bam_counts(fp, chrom_list, chrom_sizes, args.bin_size, pe_count=args.pe_count, include_duplicates=args.include_duplicates, mapq=args.mapq, norm=args.norm)
            all_vals = np.concatenate([vals[c] for c in chrom_list]) if chrom_list else np.array([0.0])
            sc = scale_values(all_vals, method=args.scale_method)
            tracks.append({"name": name, "color": color, "values_by_chrom": vals, "scale": sc, "clip_val": make_clip_val(sc), "style": args.rna_style})

    if not tracks:
        raise SystemExit("ERROR: No tracks generated.")

    gene_mid_by_chrom = None
    if args.gene_bed:
        gb = os.path.expanduser(args.gene_bed)
        if not os.path.exists(gb):
            raise SystemExit(f"ERROR: gene bed not found: {gb}")
        gene_mid_by_chrom = read_gene_midpoints(gb, chrom_list)

    out_prefix = os.path.expanduser(args.out_prefix)
    out_dir = os.path.dirname(out_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    out_fig = f"{out_prefix}.circos.{args.fig_format}"

    plot_circos(
        tracks=tracks,
        chrom_list=chrom_list,
        chrom_start=chrom_start,
        chrom_nbins=chrom_nbins,
        total_units=total_units,
        out_fig=out_fig,
        bin_size=args.bin_size,
        gene_mid_by_chrom=gene_mid_by_chrom,
        gene_tick_color=args.gene_tick_color,
        gene_tick_len=args.gene_tick_len,
        gene_tick_lw=args.gene_tick_lw,
        gene_tick_stride=args.gene_tick_stride,
        fig_size=args.fig_size,
        dpi=args.dpi,
        ideogram_radius=args.ideogram_radius,
        ideogram_width=args.ideogram_width,
        chr_label_radius=args.chr_label_radius,
        base_radius=args.base_radius,
        track_height=args.track_height,
        track_gap=args.track_gap,
        rotate_deg=args.rotate_deg,
        outline_lw=args.outline_lw,
        seg_lw=args.seg_lw,
        fill_alpha=args.fill_alpha,
        legend=(not args.no_legend),
        legend_fontsize=args.legend_fontsize,
        ytick_angle_deg=args.ytick_angle_deg,
    )

    print(f"[OK] Output figure: {out_fig}")


if __name__ == "__main__":
    main()
