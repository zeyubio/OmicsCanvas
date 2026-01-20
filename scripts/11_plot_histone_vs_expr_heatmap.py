#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
omicscanvas_histone_vs_expr_heatmap_EN.py
======================================

Purpose
-------
Plot ChIP-seq/ATAC-seq (or other track) matrix heatmaps binned by RNA expression.

This script is designed to work **directly with the matrices produced by**
`omicscanvas_bam_to_gene_matrices.py` (Step2). It assumes the matrix columns
are already ordered along the gene axis (no `po` shifting/concatenation needed).

Supported gene regions (gene_type)
----------------------------------
- ``TSS``: windows around transcription start site
- ``gene``: upstream + scaled gene body + downstream
- ``TES``: windows around transcription end site

File naming convention (standard)
--------------------------------
By default, the script reads the **new standardized suffixes**:

- ``*_tss_matrix.tsv``
- ``*_gene_profile_matrix.tsv``
- ``*_tes_matrix.tsv``

If you still have legacy outputs, override with ``--suffix-tss/--suffix-gene/--suffix-tes``.

Tracks / replicates syntax
--------------------------
Use ``--tracks`` to specify marks. Comma-separated items.

- ``NAME``: display name = NAME; file id = ``<file-prefix>NAME``
- ``NAME:FILEID``: display name = NAME; file id = ``<file-prefix>FILEID``
- Replicates: ``NAME:REP1+REP2`` (replicates are averaged)

Examples:
  --tracks "ATAC,H3K4me1,H3K4me3" --file-prefix "B_"
  --tracks "ATAC:ATAC_rep1+ATAC_rep2,H3K4me3" --file-prefix "B_"

Expression binning
------------------
- genes with expression <= ``--zero-threshold`` are split into ``--none-bins`` groups
- genes with expression >  ``--zero-threshold`` are split into ``--exp-bins`` groups
Each group contributes one heatmap row (mean/median across genes).

Compatibility
-------------
Python 3.8+.
"""

import os
import argparse
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# -------------------------
# Helpers
# -------------------------

def ensure_dir(d: str) -> str:
    return d if d.endswith("/") else (d + "/")


def split_csv(s: str) -> List[str]:
    return [x.strip() for x in s.split(",") if x.strip()]


def parse_int_list(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_float_pair(s: str) -> Tuple[float, float]:
    a, b = [x.strip() for x in s.split(",")]
    return float(a), float(b)


def is_number_like(x) -> bool:
    try:
        float(str(x))
        return True
    except Exception:
        return False


def sort_columns_numeric(cols: List[str]) -> List[str]:
    """Sort bin columns numerically if possible; otherwise keep original order."""
    if all(is_number_like(c) for c in cols):
        return sorted(cols, key=lambda z: float(str(z)))
    return cols


def safe_name(s: str) -> str:
    """Convert a label to a safe file-name token."""
    return "".join([c if (c.isalnum() or c in "._-") else "_" for c in s])


# -------------------------
# vmin / vmax preset
# -------------------------

def default_cmap(scale_mode: str) -> str:
    if scale_mode == "ratio":
        return "viridis"
    return "RdBu_r"


def compute_vmin_vmax(
    heat: np.ndarray,
    scale_mode: str,
    user_range: Optional[float],
    user_vmin: Optional[float],
    user_vmax: Optional[float],
    q_low: float,
    q_high: float,
) -> Tuple[float, float]:
    """Infer (vmin, vmax) given data + user hints.

    Priority:
      1) user_vmin and user_vmax provided -> use directly
      2) user_range provided -> ratio/quantile: [0, range]; diverging: [-range, +range]
      3) quantile -> use [q_low, q_high] quantiles
      4) diverging -> symmetric ±abs(99% quantile)
      5) ratio -> vmin=0, vmax=99% quantile
    """
    data = heat[np.isfinite(heat)]
    if data.size == 0:
        return 0.0, 1.0

    if (user_vmin is not None) and (user_vmax is not None):
        vmin = float(user_vmin)
        vmax = float(user_vmax)
        if vmax <= vmin:
            vmax = vmin + 1e-6
        return vmin, vmax

    if user_range is not None:
        r = float(abs(user_range))
        if scale_mode == "diverging":
            vmin = -r if user_vmin is None else float(user_vmin)
            vmax = r if user_vmax is None else float(user_vmax)
        else:
            vmin = 0.0 if user_vmin is None else float(user_vmin)
            vmax = r if user_vmax is None else float(user_vmax)
        if vmax <= vmin:
            vmax = vmin + 1e-6
        return vmin, vmax

    if scale_mode == "quantile":
        lo = float(np.nanquantile(data, q_low))
        hi = float(np.nanquantile(data, q_high))
        vmin = lo if user_vmin is None else float(user_vmin)
        vmax = hi if user_vmax is None else float(user_vmax)
        if vmax <= vmin:
            vmax = vmin + 1e-6
        return vmin, vmax

    if scale_mode == "diverging":
        peak = float(np.nanquantile(np.abs(data), 0.99))
        vmin = -peak if user_vmin is None else float(user_vmin)
        vmax = peak if user_vmax is None else float(user_vmax)
        if vmax <= vmin:
            vmax = vmin + 1e-6
        return vmin, vmax

    # ratio
    vmin = 0.0 if user_vmin is None else float(user_vmin)
    vmax = float(np.nanquantile(data, 0.99)) if user_vmax is None else float(user_vmax)
    if vmax <= vmin:
        vmax = vmin + 1e-6
    return vmin, vmax


# -------------------------
# Gene ID matching
# -------------------------

def normalize_gene_id(g: str, mode: str) -> str:
    g = str(g)
    if mode == "exact":
        return g
    if mode == "strip_dot":
        return g.split(".")[0]
    raise ValueError(f"Unknown gene-id-mode: {mode}")


def collapse_by_gene(df: pd.DataFrame, mode: str, how: str) -> pd.DataFrame:
    """After strip_dot, multiple isoforms may map to one gene; collapse them."""
    if mode == "exact":
        return df

    df2 = df.copy()
    df2.index = [normalize_gene_id(x, mode) for x in df2.index]

    if how == "first":
        return df2[~df2.index.duplicated(keep="first")]
    if how == "max":
        return df2.groupby(df2.index).max()
    return df2.groupby(df2.index).mean()


# -------------------------
# Read expression
# -------------------------

def load_expression(
    expr_path: str,
    expr_cols: List[int],
    expr_name: str,
    gene_id_mode: str,
    collapse: str,
) -> pd.DataFrame:
    if not os.path.exists(expr_path):
        raise FileNotFoundError(f"Not found: {expr_path}")

    df = pd.read_csv(expr_path, sep="\t", index_col=0)
    if df.shape[1] == 0:
        raise ValueError("Expression table has no columns.")

    max_col = df.shape[1] - 1
    for c in expr_cols:
        if c < 0 or c > max_col:
            raise ValueError(f"expr-col {c} out of range (0..{max_col})")

    v = pd.to_numeric(df.iloc[:, expr_cols].mean(axis=1), errors="coerce").fillna(0.0)
    out = pd.DataFrame({expr_name: v})
    out.index = out.index.astype(str)
    out = collapse_by_gene(out, gene_id_mode, collapse)
    return out


def bin_genes_by_expression(
    expr_df: pd.DataFrame,
    expr_name: str,
    none_bins: int,
    exp_bins: int,
    zero_threshold: float,
    ascending: bool,
) -> List[List[str]]:
    s = expr_df[expr_name].sort_values(ascending=ascending)
    none_part = s[s <= zero_threshold]
    exp_part = s[s > zero_threshold]

    groups: List[List[str]] = []

    if none_bins > 0:
        splits = (
            np.array_split(none_part.index.to_numpy(), none_bins)
            if len(none_part) > 0
            else [np.array([]) for _ in range(none_bins)]
        )
        groups.extend([arr.tolist() for arr in splits])

    if exp_bins > 0:
        splits = (
            np.array_split(exp_part.index.to_numpy(), exp_bins)
            if len(exp_part) > 0
            else [np.array([]) for _ in range(exp_bins)]
        )
        groups.extend([arr.tolist() for arr in splits])

    return groups


# -------------------------
# Read matrices
# -------------------------

def read_matrix(
    path: str,
    gene_id_mode: str,
    collapse: str,
    isoform_suffix: Optional[str],
) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")

    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)

    if isoform_suffix:
        df = df.loc[[g for g in df.index if g.endswith(isoform_suffix)]].copy()

    num = df.apply(pd.to_numeric, errors="coerce")
    num = num.loc[:, num.notna().any(axis=0)]

    cols = sort_columns_numeric(list(num.columns))
    num = num.loc[:, cols]

    num = collapse_by_gene(num, gene_id_mode, collapse)
    return num


def average_replicates(
    file_ids: List[str],
    chip_dir: str,
    suffix: str,
    gene_id_mode: str,
    collapse: str,
    isoform_suffix: Optional[str],
) -> pd.DataFrame:
    base = ensure_dir(chip_dir)
    mats: List[pd.DataFrame] = []
    for fid in file_ids:
        path = os.path.join(base, fid + suffix)
        mats.append(read_matrix(path, gene_id_mode, collapse, isoform_suffix))

    if len(mats) == 1:
        return mats[0]

    common_cols = set(mats[0].columns)
    for m in mats[1:]:
        common_cols &= set(m.columns)
    if len(common_cols) == 0:
        raise ValueError("No common bins(columns) across replicates.")
    common_cols = sort_columns_numeric(list(common_cols))

    all_genes = sorted(set().union(*[set(m.index) for m in mats]))

    stack = np.stack(
        [
            m.reindex(index=all_genes, columns=common_cols).to_numpy(dtype=float)
            for m in mats
        ],
        axis=0,
    )

    avg = np.nanmean(stack, axis=0)
    return pd.DataFrame(avg, index=all_genes, columns=common_cols)


# -------------------------
# Tracks parsing
# -------------------------

def parse_tracks(tracks_str: str, file_prefix: str) -> List[Tuple[str, List[str]]]:
    """Parse --tracks specification."""
    tracks: List[Tuple[str, List[str]]] = []
    file_prefix = file_prefix or ""

    for item in split_csv(tracks_str):
        if ":" in item:
            display, fid_part = item.split(":", 1)
            display = display.strip()
            fid_part = fid_part.strip()
        else:
            display = item.strip()
            fid_part = item.strip()

        rep_ids = [x.strip() for x in fid_part.split("+") if x.strip()]

        rep_ids2: List[str] = []
        for rid in rep_ids:
            if file_prefix and (not rid.startswith(file_prefix)):
                rep_ids2.append(file_prefix + rid)
            else:
                rep_ids2.append(rid)

        if not display:
            raise ValueError(f"Invalid track item: {item}")
        if len(rep_ids2) == 0:
            raise ValueError(f"No file id parsed for track: {item}")

        tracks.append((display, rep_ids2))

    return tracks


# -------------------------
# Heat matrix
# -------------------------

def build_heat(avg_mat: pd.DataFrame, gene_groups: List[List[str]], stat: str) -> np.ndarray:
    cols = list(avg_mat.columns)
    n_bins = len(cols)

    out_rows: List[np.ndarray] = []
    for genes in gene_groups:
        if len(genes) == 0:
            out_rows.append(np.full(n_bins, np.nan))
            continue

        g = [x for x in genes if x in avg_mat.index]
        if len(g) == 0:
            out_rows.append(np.full(n_bins, np.nan))
            continue

        sub = avg_mat.loc[g, cols].to_numpy(dtype=float)
        if stat == "median":
            out_rows.append(np.nanmedian(sub, axis=0))
        else:
            out_rows.append(np.nanmean(sub, axis=0))

    return np.vstack(out_rows)


# -------------------------
# Axis templates
# -------------------------

def axis_template(
    gene_type: str,
    distance: int,
    bins_start: int,
    bins_end: int,
    bins_gene_st: int,
    bins_gene_body: int,
    bins_gene_en: int,
) -> Tuple[int, List[int], List[str], List[int]]:
    """Return: bins_total, xtick_positions, xtick_labels, vertical_line_positions."""
    if gene_type == "TSS":
        bins_total = bins_start
        mid = bins_start // 2 - 1
        xtick_pos = [0, mid, bins_total - 1]
        xtick_lab = [f"-{distance}bp", "TSS", f"+{distance}bp"]
        xline_pos = [mid]
        return bins_total, xtick_pos, xtick_lab, xline_pos

    if gene_type == "TES":
        bins_total = bins_end
        mid = bins_end // 2 - 1
        xtick_pos = [0, mid, bins_total - 1]
        xtick_lab = [f"-{distance}bp", "TES", f"+{distance}bp"]
        xline_pos = [mid]
        return bins_total, xtick_pos, xtick_lab, xline_pos

    bins_total = bins_gene_st + bins_gene_body + bins_gene_en
    tss = bins_gene_st - 1
    tes = bins_gene_st + bins_gene_body - 1
    xtick_pos = [0, tss, tes, bins_total - 1]
    xtick_lab = [f"-{distance}bp", "TSS", "TES", f"+{distance}bp"]
    xline_pos = [tss, tes]
    return bins_total, xtick_pos, xtick_lab, xline_pos


# -------------------------
# Plot
# -------------------------

def plot_heatmap(
    heat: np.ndarray,
    out_file: str,
    xtick_pos: List[int],
    xtick_lab: List[str],
    xline_pos: List[int],
    cmap: str,
    vmin: float,
    vmax: float,
    ytick_step: Optional[int],
    border: bool,
    title: Optional[str],
    dpi: int,
    out_format: str,
):
    nrows, ncols = heat.shape

    fig_h = max(6, nrows / 10.0)
    fig_w = 8
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = plt.gca()

    sns.heatmap(heat, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, cbar=True)

    ax.set_xticks([x + 0.5 for x in xtick_pos])
    ax.set_xticklabels(xtick_lab, rotation=0)

    if ytick_step is None:
        ytick_step = max(1, nrows // 10)
    yticks = list(range(0, nrows, ytick_step))
    ax.set_yticks([y + 0.5 for y in yticks])
    ax.set_yticklabels([str(y) for y in yticks], rotation=0)

    for xp in xline_pos:
        ax.axvline(xp, linestyle="--", color="black", linewidth=1)

    if border:
        ax.plot([0, ncols], [0, 0], color="black", lw=1)
        ax.plot([0, ncols], [nrows, nrows], color="black", lw=1)
        ax.plot([0, 0], [0, nrows], color="black", lw=1)
        ax.plot([ncols, ncols], [0, nrows], color="black", lw=1)

    if title:
        ax.set_title(title)

    plt.savefig(out_file, bbox_inches="tight", dpi=dpi, format=out_format)
    plt.close(fig)


# -------------------------
# Main
# -------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Build heatmaps of track matrices (ChIP/ATAC) after binning genes by expression. "
            "Input matrices are generated by omicscanvas_bam_to_gene_matrices.py."
        ),
    )

    # matrix inputs
    ap.add_argument(
        "--matrix-dir",
        required=True,
        help="Directory containing per-track matrices (one file per mark per region).",
    )
    ap.add_argument(
        "--tracks",
        required=True,
        help=(
            "Track list (comma-separated). Each item: NAME or NAME:FILEID or NAME:REP1+REP2. "
            "FILEID(s) will be prefixed by --file-prefix if they do not already start with it."
        ),
    )
    ap.add_argument(
        "--file-prefix",
        default="",
        help="Common prefix to add to each file id, e.g. 'B_'.",
    )

    # gene types & suffix
    ap.add_argument(
        "--gene-types",
        default="gene",
        help="Which regions to plot: TSS,gene,TES (comma-separated).",
    )
    ap.add_argument(
        "--suffix-tss",
        default="_tss_matrix.tsv",
        help="Suffix for TSS matrix files (standard: _tss_matrix.tsv).",
    )
    ap.add_argument(
        "--suffix-gene",
        default="_gene_profile_matrix.tsv",
        help="Suffix for gene-profile matrix files (standard: _gene_profile_matrix.tsv).",
    )
    ap.add_argument(
        "--suffix-tes",
        default="_tes_matrix.tsv",
        help="Suffix for TES matrix files (standard: _tes_matrix.tsv).",
    )

    # bins and distance (ticks/lines + sanity check)
    ap.add_argument("--distance", type=int, default=2000, help="Distance shown on x-axis labels (bp).")
    ap.add_argument("--bins-start", type=int, default=100, help="Total bins for TSS region.")
    ap.add_argument("--bins-end", type=int, default=100, help="Total bins for TES region.")
    ap.add_argument("--bins-gene-st", type=int, default=50, help="Gene profile: upstream bins.")
    ap.add_argument("--bins-gene-body", type=int, default=100, help="Gene profile: scaled gene-body bins.")
    ap.add_argument("--bins-gene-en", type=int, default=50, help="Gene profile: downstream bins.")

    # gene id matching
    ap.add_argument(
        "--gene-id-mode",
        choices=["exact", "strip_dot"],
        default="exact",
        help=(
            "How to match gene IDs between expression and matrices. "
            "exact: keep full IDs; strip_dot: remove isoform suffix after '.' (e.g. AT1G01010.1 -> AT1G01010)."
        ),
    )
    ap.add_argument(
        "--collapse-isoform",
        choices=["mean", "max", "first"],
        default="mean",
        help="If gene-id-mode=strip_dot, how to collapse multiple isoforms per gene.",
    )
    ap.add_argument(
        "--isoform-suffix",
        default=None,
        help="Optional: keep only matrix rows ending with this suffix (e.g. '.1').",
    )

    # expression
    ap.add_argument("--expr", required=True, help="Expression table (TSV). First column is gene ID (index).")
    ap.add_argument(
        "--expr-cols",
        default="0,1,2",
        help="0-based column indices (after the index column) used to compute mean expression.",
    )
    ap.add_argument("--expr-name", default="EXP", help="Internal expression column name.")

    # binning
    ap.add_argument("--none-bins", type=int, default=10, help="Number of bins for low/zero-expression genes.")
    ap.add_argument("--exp-bins", type=int, default=90, help="Number of bins for expressed genes.")
    ap.add_argument(
        "--zero-threshold",
        type=float,
        default=0.0,
        help="Expression <= this value is treated as 'none/low' group.",
    )
    ap.add_argument(
        "--descending",
        action="store_true",
        help="Sort genes by expression high->low (default is low->high).",
    )

    # statistic within each expression bin
    ap.add_argument(
        "--stat",
        choices=["mean", "median"],
        default="mean",
        help="Statistic used within each expression bin to form one heatmap row.",
    )

    # color scaling
    ap.add_argument(
        "--scale-mode",
        choices=["quantile", "diverging", "ratio"],
        default="quantile",
        help=(
            "Auto color scaling mode. quantile: use --quantiles; diverging: symmetric; ratio: non-negative (0..max)."
        ),
    )
    ap.add_argument(
        "--quantiles",
        default="0.01,0.99",
        help="Quantiles used when scale-mode=quantile, e.g. 0.05,0.95.",
    )
    ap.add_argument(
        "--range",
        type=float,
        default=None,
        help="Provide one number to set range: diverging -> ±range; otherwise -> [0, range].",
    )
    ap.add_argument("--vmin", type=float, default=None, help="Manual vmin (overrides auto).")
    ap.add_argument("--vmax", type=float, default=None, help="Manual vmax (overrides auto).")
    ap.add_argument("--cmap", default=None, help="Colormap name (e.g. RdBu_r, viridis).")

    # plot style
    ap.add_argument("--ytick-step", type=int, default=None, help="Y tick step (default auto).")
    ap.add_argument("--no-border", action="store_true", help="Disable outer border box.")
    ap.add_argument("--title", default=None, help="Optional title for each heatmap.")
    ap.add_argument("--dpi", type=int, default=300, help="Output DPI (for PNG/PDF rasterization).")
    ap.add_argument(
        "--out-format",
        choices=["pdf", "png"],
        default="pdf",
        help="Output format.",
    )

    # output
    ap.add_argument(
        "--out-prefix",
        default="out",
        help="Output prefix. Files will be named <out-prefix>_<track>_<type>_histone_vs_expr_heatmap.<ext>",
    )
    ap.add_argument(
        "--outdir",
        default=".",
        help="Output directory.",
    )

    args = ap.parse_args()

    gene_types = [x.strip() for x in split_csv(args.gene_types) if x.strip()]
    expr_cols = parse_int_list(args.expr_cols)

    tracks = parse_tracks(args.tracks, args.file_prefix)

    expr_df = load_expression(
        expr_path=args.expr,
        expr_cols=expr_cols,
        expr_name=args.expr_name,
        gene_id_mode=args.gene_id_mode,
        collapse=args.collapse_isoform,
    )

    q_low, q_high = parse_float_pair(args.quantiles)
    cmap = args.cmap if args.cmap is not None else default_cmap(args.scale_mode)

    os.makedirs(args.outdir, exist_ok=True)

    suffix_map = {"TSS": args.suffix_tss, "gene": args.suffix_gene, "TES": args.suffix_tes}

    for gtype in gene_types:
        if gtype not in ("TSS", "gene", "TES"):
            raise ValueError(f"Invalid gene-type: {gtype}. Use TSS,gene,TES")

        suffix = suffix_map[gtype]

        expected_bins, xtick_pos, xtick_lab, xline_pos = axis_template(
            gtype,
            distance=args.distance,
            bins_start=args.bins_start,
            bins_end=args.bins_end,
            bins_gene_st=args.bins_gene_st,
            bins_gene_body=args.bins_gene_body,
            bins_gene_en=args.bins_gene_en,
        )

        for display_name, file_ids in tracks:
            mat = average_replicates(
                file_ids=file_ids,
                chip_dir=args.matrix_dir,
                suffix=suffix,
                gene_id_mode=args.gene_id_mode,
                collapse=args.collapse_isoform,
                isoform_suffix=args.isoform_suffix,
            )

            n_bins = mat.shape[1]
            if n_bins != expected_bins:
                raise ValueError(
                    f"[{gtype}|{display_name}] bins mismatch: matrix has {n_bins} bins, "
                    f"but expected {expected_bins} from --bins-* settings."
                )

            overlap = sorted(set(expr_df.index) & set(mat.index))
            if len(overlap) == 0:
                raise ValueError(
                    f"No overlap genes between expression and matrix: type={gtype} track={display_name}"
                )

            expr_overlap = expr_df.loc[overlap]
            groups = bin_genes_by_expression(
                expr_overlap,
                expr_name=args.expr_name,
                none_bins=args.none_bins,
                exp_bins=args.exp_bins,
                zero_threshold=args.zero_threshold,
                ascending=(not args.descending),
            )

            heat = build_heat(mat, groups, stat=args.stat)

            vmin, vmax = compute_vmin_vmax(
                heat=heat,
                scale_mode=args.scale_mode,
                user_range=args.range,
                user_vmin=args.vmin,
                user_vmax=args.vmax,
                q_low=q_low,
                q_high=q_high,
            )

            mark_safe = safe_name(display_name)
            out_ext = args.out_format
            out_file = os.path.join(
                args.outdir,
                f"{args.out_prefix}_{mark_safe}_{gtype}_histone_vs_expr_heatmap.{out_ext}",
            )

            plot_heatmap(
                heat=heat,
                out_file=out_file,
                xtick_pos=xtick_pos,
                xtick_lab=xtick_lab,
                xline_pos=xline_pos,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                ytick_step=args.ytick_step,
                border=(not args.no_border),
                title=args.title,
                dpi=args.dpi,
                out_format=args.out_format,
            )

            print(
                f"[OK] {out_file} | type={gtype} track={display_name} files={'+'.join(file_ids)} "
                f"cmap={cmap} vmin={vmin:.6g} vmax={vmax:.6g} overlap={len(overlap)}"
            )


if __name__ == "__main__":
    main()
