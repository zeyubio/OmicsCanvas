#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
17_methylation_bin_correlation_triangle.py

Bin methylation CX files (chrom, pos, meth_count, depth) by fixed bp bins,
compute per-bin methylation ratio (sum(meth)/sum(depth)),
intersect bins across samples, compute sample-sample correlation,
and plot a *triangular* correlation heatmap (diamond cells) like the example figure:
- draw ONLY one triangle (default: lower triangle, i>j)
- remove diagonal (self-self = 1) by not plotting it
- orient the longest edge (base) as X axis (horizontal)

Input format (tab-delimited, no header):
Chr01   15072   1   1
Chr01   15073   0   1
Chr01   15347   8   8
Chr01   15348   4   4
Columns:
1) chrom
2) pos (1-based)
3) methylated count
4) depth

Outputs:
- <out-prefix>.<group>.corr.tsv : correlation matrix (samples x samples)
- <out-prefix>.corr.triangle.<pdf|png> : triangular correlation figure (single or multi-panel)
- optional <out-prefix>.<group>.bin_matrix.tsv : per-bin methylation ratio (bins x samples)

Multi-panel usage:
Provide multiple file groups separated by "|" to plot stacked panels (e.g., CG|CHG|CHH):
  --files "S2_CG.cx,S4_CG.cx,S6_CG.cx|S2_CHG.cx,S4_CHG.cx,S6_CHG.cx|S2_CHH.cx,S4_CHH.cx,S6_CHH.cx"
  --group-names "CG|CHG|CHH"
  --names "S2,S4,S6"   (optional; reuse same sample labels for each group)

Author: generated from user's request (triangle-only correlation plot).
"""

import os
import sys
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42   # TrueType
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['svg.fonttype'] = 'none'  # SVG 中保留文字
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 14



def compute_bin_ratio(
    filepath: str,
    bin_size: int = 10000,
    min_bin_depth: int = 1,
    chunksize: int = 2_000_000,
) -> pd.Series:
    """Aggregate CX-like methylation file into (chrom, bin)->ratio."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")

    meth_sum = defaultdict(int)
    depth_sum = defaultdict(int)

    colnames = ["chrom", "pos", "meth", "depth"]
    reader = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=colnames,
        usecols=[0, 1, 2, 3],
        dtype={0: "string", 1: "int64", 2: "int64", 3: "int64"},
        comment="#",
        chunksize=chunksize,
        engine="c",
    )

    for chunk in reader:
        chunk = chunk.dropna()
        chunk = chunk[chunk["depth"] > 0]
        if chunk.empty:
            continue

        chunk["bin"] = (chunk["pos"] - 1) // bin_size

        grp = chunk.groupby(["chrom", "bin"], sort=False)[["meth", "depth"]].sum()

        for (chrom, b), row in grp.iterrows():
            key = (str(chrom), int(b))
            meth_sum[key] += int(row["meth"])
            depth_sum[key] += int(row["depth"])

    if not meth_sum:
        return pd.Series(dtype="float64")

    idx = []
    vals = []
    for key in meth_sum.keys():
        d = depth_sum[key]
        if d >= min_bin_depth:
            idx.append(key)
            vals.append(meth_sum[key] / d)

    s = pd.Series(
        vals,
        index=pd.MultiIndex.from_tuples(idx, names=["chrom", "bin"]),
        dtype="float64",
    ).sort_index()
    s.name = os.path.basename(filepath)
    return s


def choose_vrange(
    corr: pd.DataFrame,
    mode: str,
    vpercentile: str,
    vmin_arg,
    vmax_arg,
) -> tuple[float, float]:
    """Choose vmin/vmax for heatmap based on off-diagonal values."""
    n = corr.shape[0]
    vals = corr.values.astype(float, copy=True)

    mask = ~np.eye(n, dtype=bool)
    offdiag = vals[mask]
    offdiag = offdiag[np.isfinite(offdiag)]

    if (vmin_arg is not None) and (vmax_arg is not None):
        vmin, vmax = float(vmin_arg), float(vmax_arg)
        if np.isfinite(vmin) and np.isfinite(vmax) and vmin < vmax:
            return vmin, vmax
        return -1.0, 1.0

    if offdiag.size == 0:
        return -1.0, 1.0

    if mode == "fixed":
        vmin, vmax = -1.0, 1.0
    elif mode == "data":
        vmin, vmax = float(np.min(offdiag)), float(np.max(offdiag))
    else:  # percentile
        try:
            p_low, p_high = [float(x) for x in vpercentile.split(",")]
        except Exception:
            p_low, p_high = 5.0, 95.0
        vmin, vmax = np.nanpercentile(offdiag, [p_low, p_high])
        vmin, vmax = float(vmin), float(vmax)

    if not (np.isfinite(vmin) and np.isfinite(vmax) and vmin < vmax):
        return -1.0, 1.0
    return vmin, vmax


def plot_corr_triangle(
    ax,
    corr: pd.DataFrame,
    labels: list[str],
    cmap: str = "viridis",
    vmin: float = 0.0,
    vmax: float = 1.0,
    triangle: str = "lower",   # "lower" (i>j) or "upper" (i<j)
    show_labels: bool = True,
    label_fontsize: int = 10,
    edgecolor: str = "white",
    linewidth: float = 0.8,
    annotate: bool = False,
    annot_fmt: str = ".2f",
    annot_fontsize: int = 9,
):
    """
    Draw a diamond-cell triangular heatmap on `ax`.
    Coordinate system:
      For lower triangle (i>j), place cell center at:
        cx = i + j
        cy = i - j
      This makes the longest edge (base) horizontal (X axis), like the example.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize

    n = corr.shape[0]
    if n < 2:
        ax.text(0.5, 0.5, "Need >=2 samples", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return None

    # prepare polygons
    polys = []
    vals = []
    centers = []

    for i in range(n):
        for j in range(n):
            if triangle == "lower":
                # keep diagonal: plot only i>=j
                if i < j:
                    continue
            else:  # upper
                # keep diagonal: plot only i<=j
                if i > j:
                    continue

            val = float(corr.iat[i, j])
            if not np.isfinite(val):
                continue

            cx = i + j
            cy = (i - j) if triangle == "lower" else (j - i)

            poly = [(cx, cy - 1), (cx + 1, cy), (cx, cy + 1), (cx - 1, cy)]
            polys.append(poly)
            vals.append(val)
            centers.append((cx, cy))

    if not polys:
        ax.text(0.5, 0.5, "No off-diagonal cells", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return None

    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    coll = PolyCollection(
        polys,
        array=np.array(vals, dtype=float),
        cmap=plt.get_cmap(cmap),
        norm=norm,
        edgecolors=edgecolor,
        linewidths=linewidth,
    )
    ax.add_collection(coll)
    ax.set_aspect("equal")

    # limits
    xmax = 2 * (n - 1)
    ymax = (n - 1)
    ax.set_xlim(-1, xmax + 1)
    ax.set_ylim(-1.5, ymax + 1.5)

    # labels on two slanted edges
    if show_labels and labels and len(labels) == n:
        # left edge: x=y, bottom->top is labels[0..n-1]
        for k, lab in enumerate(labels):
            ax.text(
                k - 0.15,
                k + 0.15,
                lab,
                rotation=45,
                ha="right",
                va="bottom",
                fontsize=label_fontsize,
            )
        # right edge: x=2*(n-1)-y, bottom->top is reversed labels
        for k, lab in enumerate(labels[::-1]):
            x = 2 * (n - 1) - k
            y = k
            ax.text(
                x + 0.15,
                y + 0.15,
                lab,
                rotation=-45,
                ha="left",
                va="bottom",
                fontsize=label_fontsize,
            )

    # optional numbers inside diamonds
    if annotate:
        for (cx, cy), v in zip(centers, vals):
            ax.text(cx, cy, format(v, annot_fmt), ha="center", va="center", fontsize=annot_fontsize)

    ax.axis("off")
    return coll


def parse_groups(files_arg: str) -> list[list[str]]:
    """Parse --files: groups separated by '|' and items by ','."""
    groups = []
    for g in files_arg.split("|"):
        g = g.strip()
        if not g:
            continue
        items = [x.strip() for x in g.split(",") if x.strip()]
        if items:
            groups.append(items)
    return groups


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Compute binned methylation ratio from CX files, intersect common bins, "
            "calculate sample correlation, and plot *triangular* correlation heatmap (diamond cells)."
        ),
    )
    ap.add_argument("--indir", required=True, help="Input directory containing CX files.")
    ap.add_argument(
        "--files",
        required=True,
        help=(
            "Files for correlation. "
            "Use comma-separated list for one panel; use '|' to separate multiple panels. "
            'Example: "S2_CG.cx,S4_CG.cx,S6_CG.cx|S2_CHG.cx,S4_CHG.cx,S6_CHG.cx|S2_CHH.cx,S4_CHH.cx,S6_CHH.cx"'
        ),
    )
    ap.add_argument(
        "--group-names",
        default=None,
        help='Optional panel names separated by "|", e.g. "CG|CHG|CHH".',
    )
    ap.add_argument(
        "--names",
        default=None,
        help=(
            "Optional sample labels (comma-separated) to use on both edges. "
            "If not provided, uses file basenames."
        ),
    )
    ap.add_argument("--bin-size", type=int, default=10000, help="Bin size in bp.")
    ap.add_argument(
        "--method",
        choices=["pearson", "spearman", "kendall"],
        default="pearson",
        help="Correlation method.",
    )
    ap.add_argument(
        "--min-bin-depth",
        type=int,
        default=10,
        help="Filter out bins whose total depth (sum depth in bin) < this threshold.",
    )
    ap.add_argument(
        "--chunksize",
        type=int,
        default=2_000_000,
        help="Read file in chunks to reduce memory usage.",
    )
    ap.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix, e.g. results/CG_CHG_CHH_10kb",
    )
    ap.add_argument(
        "--save-bin-matrix",
        action="store_true",
        help="Also save per-bin methylation ratio matrix (bins x samples) for each panel.",
    )

    ap.add_argument("--fig-format", choices=["pdf", "png"], default="pdf", help="Figure format.")
    ap.add_argument("--dpi", type=int, default=300, help="DPI for PNG (ignored for PDF).")

    ap.add_argument("--cmap", default="viridis", help='Colormap name (e.g. "viridis", "RdBu_r").')
    ap.add_argument(
        "--vrange",
        choices=["fixed", "data", "percentile"],
        default="percentile",
        help="How to choose vmin/vmax based on off-diagonal correlations.",
    )
    ap.add_argument("--vmin", type=float, default=None, help="Manual vmin (effective only if --vmax also set).")
    ap.add_argument("--vmax", type=float, default=None, help="Manual vmax (effective only if --vmin also set).")
    ap.add_argument("--vpercentile", default="5,95", help='Percentiles for vrange=percentile, e.g. "5,95".')

    ap.add_argument("--triangle", choices=["lower", "upper"], default="lower", help="Which triangle to plot.")
    ap.add_argument("--annot", action="store_true", help="Annotate each diamond with correlation value.")

    args = ap.parse_args()

    # non-interactive backend
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    indir = args.indir
    groups = parse_groups(args.files)
    if not groups:
        raise ValueError("No valid files parsed from --files.")
    n_panels = len(groups)

    group_names = None
    if args.group_names:
        group_names = [x.strip() for x in args.group_names.split("|")]
        if len(group_names) != n_panels:
            print("[WARN] --group-names count != number of panels; ignoring names.", file=sys.stderr)
            group_names = None
    if group_names is None:
        group_names = [f"group{i+1}" for i in range(n_panels)]

    # sample labels
    user_labels = None
    if args.names:
        user_labels = [x.strip() for x in args.names.split(",") if x.strip()]

    # output dir
    out_prefix = args.out_prefix
    out_dir = os.path.dirname(out_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # layout: each panel with its own colorbar beneath
    height_per_panel = 2.4
    fig_h = max(2.5, height_per_panel * n_panels)
    fig_w = 4.2
    fig = plt.figure(figsize=(fig_w, fig_h), constrained_layout=True)

    # gridspec with (panel + colorbar) rows repeated
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(nrows=n_panels * 2, ncols=1, figure=fig, height_ratios=[1.0, 0.08] * n_panels)

    for p, (gname, files) in enumerate(zip(group_names, groups)):
        if len(files) < 2:
            raise ValueError(f"Panel {gname}: need at least 2 files.")
        # compute bin ratios
        series_list = []
        sample_names = []
        for fn in files:
            fp = os.path.join(indir, fn)
            s = compute_bin_ratio(fp, bin_size=args.bin_size, min_bin_depth=args.min_bin_depth, chunksize=args.chunksize)
            if s.empty:
                print(f"[WARN] No valid bins after filtering for file: {fp}", file=sys.stderr)
            series_list.append(s)
            sample_names.append(os.path.basename(fn))

        bin_matrix = pd.concat(series_list, axis=1, join="inner")
        bin_matrix.columns = sample_names

        if bin_matrix.shape[0] == 0:
            raise RuntimeError(
                f"Panel {gname}: No common (chrom,bin) across samples. "
                "Try lowering --min-bin-depth or check chromosome naming."
            )

        corr = bin_matrix.corr(method=args.method)

        # save corr / bin matrix
        corr_tsv = f"{out_prefix}.{gname}.corr.tsv"
        corr.to_csv(corr_tsv, sep="\t", float_format="%.6f")
        if args.save_bin_matrix:
            bin_tsv = f"{out_prefix}.{gname}.bin_matrix.tsv"
            bm = bin_matrix.copy()
            bm.reset_index().rename(columns={"chrom": "CHR", "bin": "BIN_ID"}).to_csv(
                bin_tsv, sep="\t", index=False, float_format="%.6f"
            )

        # choose vrange
        vmin, vmax = choose_vrange(corr, args.vrange, args.vpercentile, args.vmin, args.vmax)

        # labels for the triangle
        labels = user_labels if (user_labels is not None and len(user_labels) == corr.shape[0]) else sample_names

        ax = fig.add_subplot(gs[p * 2, 0])
        mappable = plot_corr_triangle(
            ax=ax,
            corr=corr,
            labels=labels,
            cmap=args.cmap,
            vmin=vmin,
            vmax=vmax,
            triangle=args.triangle,
            show_labels=True,
            annotate=args.annot,
        )
        ax.text(0.02, 0.95, gname, transform=ax.transAxes, ha="left", va="top", fontsize=12, fontweight="bold")

        # colorbar
        cax = fig.add_subplot(gs[p * 2 + 1, 0])
        if mappable is not None:
            cb = fig.colorbar(mappable, cax=cax, orientation="horizontal")
            cb.set_label(f"Correlation ({args.method})  vmin={vmin:.3f} vmax={vmax:.3f}", fontsize=9)
            cax.tick_params(labelsize=8)
        else:
            cax.axis("off")

        print(f"[OK] {gname}: common_bins={bin_matrix.shape[0]}  corr={corr_tsv}", file=sys.stderr)

    fig_path = f"{out_prefix}.corr.triangle.{args.fig_format}"
    if args.fig_format == "png":
        fig.savefig(fig_path, dpi=args.dpi, bbox_inches="tight")
    else:
        fig.savefig(fig_path, bbox_inches="tight")

    print(f"[OK] Figure: {fig_path}")
    print("[NOTE] If plus/minus labels look reversed in biology context, verify GRO library strandedness in IGV.")


if __name__ == "__main__":
    main()
