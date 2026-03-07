#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Whole-profile methylation plot (2D + 3D layered) with custom line colors

Input filename convention:
    <meth_dir>/<sample_id>_<CX><whole_suffix>
For example:
    gene/LSH10_CG_whole_line.txt
    gene/ST_8_CHH_whole_line.txt

File content:
    - By default, the file is treated as a single-column numeric profile (one methylation value per line)
    - If multiple columns are present, the first column is used as the methylation value

Output:
    <out_prefix>.2D.pdf

Example:
    python methylation_whole_profile_2d3d.py \
      --meth-dir gene \
      --samples LSH10,ST_8 \
      --labels  LSH10,WT \
      --contexts CG,CHG,CHH \
      --distance 2000 \
      --bins-st 50 --bins-body 100 --bins-en 50 \
      --whole-suffix _whole_line.txt \
      --mode both \
      --cmap tab10 \
      --line-colors "#FF3030,#00EE76" \
      --line-lw 2 \
      --fig-x 6 --fig-y 8 \
      --edge 0.10 --x-offset 20 --y-offset 35 \
      --out-prefix 1gene_whole_profile_methylation

Notes:
- The bin counts are used only for x-axis ticks (TSS/TES) and vertical guide lines, and must match the binning used to generate the whole-profile files.
- The legacy distance labels were shown as 'left +distance / right -distance', and this script keeps that behavior by default.
  To use the conventional labeling 'left -distance / right +distance', add --standard-distance-labels.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42   # TrueType
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['font.size'] = 14
from matplotlib.patches import ConnectionPatch



# -------------------------
# Basic helpers
# -------------------------

def _split_csv(s: str):
    """Convert 'A,B,C' to ['A', 'B', 'C']; surrounding spaces are allowed."""
    return [x.strip() for x in s.split(",") if x.strip()]


def ensure_dir_slash(d: str) -> str:
    """Ensure that the directory path ends with a slash."""
    if d and (not d.endswith("/")):
        d += "/"
    return d


def load_one_profile(path: str) -> pd.DataFrame:
    """
    Read a single profile file and return a DataFrame with one column: ['methylation'].
    - Default: no header, single-column input
    - If multiple columns are present, use the first column
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[0] == 0:
        raise ValueError(f"Empty file: {path}")

    out = pd.DataFrame({"methylation": pd.to_numeric(df.iloc[:, 0], errors="coerce")})
    if out["methylation"].isna().all():
        raise ValueError(f"All values are NaN after parsing: {path}")
    return out


def load_profiles(meth_dir: str,
                  samples: list,
                  labels: list,
                  contexts: list,
                  whole_suffix: str) -> dict:
    """
    Return structure:
        profiles[CX][label] = DataFrame(['methylation'])
    """
    if len(samples) != len(labels):
        raise ValueError(f"The number of --samples ({len(samples)}) must equal the number of --labels ({len(labels)}).")

    meth_dir = ensure_dir_slash(meth_dir)

    profiles = {}
    for cx in contexts:
        profiles[cx] = {}
        for sid, lab in zip(samples, labels):
            path = f"{meth_dir}{sid}_{cx}{whole_suffix}"
            profiles[cx][lab] = load_one_profile(path)

    return profiles


def parse_line_colors(color_str: str, n_lines: int):
    """
    Parse user-defined line colors in sample/label order.

      - Named colors: red, blue, grey
      - Hex colors: #FF3030

    Rules:
      - If one color is provided, use it for all samples
      - If fewer colors than samples are provided, colors are cycled
      - If more colors than samples are provided, only the first n_lines are used
    """
    if color_str is None:
        return None

    cols = _split_csv(color_str)
    if len(cols) == 0:
        return None

    for c in cols:
        if not mcolors.is_color_like(c):
            raise ValueError(
                f"Invalid color: {c}. Use matplotlib-compatible colors, e.g. 'red' or '#FF3030'."
            )

    if len(cols) == 1:
        return cols * n_lines
    if len(cols) < n_lines:
        cols = (cols * (n_lines // len(cols) + 1))[:n_lines]
    else:
        cols = cols[:n_lines]
    return cols


def parse_ylim(ylim_str: str, contexts: list):
    """
    Parse --ylim.
    Supported formats:
      1) Single pair: "ymin,ymax"      -> shared by all contexts
      2) Multiple pairs: "ymin,ymax;ymin,ymax;..." -> the number of pairs must equal the number of contexts
    Return a dict: ylim_dict[cx] = (ymin, ymax)
    """
    if ylim_str is None:
        return None

    s = ylim_str.strip()
    if not s:
        return None

    if ";" not in s:
        a, b = [x.strip() for x in s.split(",")]
        ymin, ymax = float(a), float(b)
        return {cx: (ymin, ymax) for cx in contexts}

    parts = [p.strip() for p in s.split(";") if p.strip()]
    if len(parts) != len(contexts):
        raise ValueError(f"The number of --ylim groups must equal the number of contexts ({len(contexts)}); current={len(parts)}")
    out = {}
    for cx, p in zip(contexts, parts):
        a, b = [x.strip() for x in p.split(",")]
        out[cx] = (float(a), float(b))
    return out


def make_xticks(distance: int,
                bins_st: int,
                bins_body: int,
                bins_en: int,
                standard_distance_labels: bool = False):
    """
    Define x-axis ticks and vertical guide lines for a gene-body-normalized whole profile.

    all_length = bins_st + bins_body + bins_en
    Position meanings:
        0                          : upstream boundary (distance from TSS)
        bins_st                    : TSS
        bins_st + bins_body        : TES
        all_length                 : downstream boundary (distance from TES)

    Distance labels:
        - Default: keep the legacy style, left +distance / right -distance
        - If standard_distance_labels=True: left -distance / right +distance
    """
    all_length = bins_st + bins_body + bins_en
    x_pos = [0, bins_st, bins_st + bins_body, all_length]

    if standard_distance_labels:
        x_labels = [f"-{distance}", "TSS", "TES", f"+{distance}"]
    else:
        # Keep compatibility with the legacy labeling convention
        x_labels = [f"+{distance}", "TSS", "TES", f"-{distance}"]

    x_lines = [bins_st, bins_st + bins_body]
    return all_length, x_pos, x_labels, x_lines


# -------------------------
# ConnectionPatch helper
# -------------------------

def connect_vertical(xpos, ax_upper, ax_lower):
    """
    Draw a vertical dashed connector between two axes (coordinates are transformed automatically).

    xpos: shared x coordinate (for example, the bin index of TSS or TES)
    ax_upper: upper axis
    ax_lower: lower axis

    Note: the lower bound of each axis (get_ylim()[0]) is used as the connector endpoint,
         consistent with the original implementation.
    """
    y_bottom_of_upper = ax_upper.get_ylim()[0]
    y_top_of_lower = ax_lower.get_ylim()[0]

    fig = ax_upper.figure
    conn = ConnectionPatch(
        xyA=(xpos, y_bottom_of_upper), coordsA=ax_upper.transData,
        xyB=(xpos, y_top_of_lower), coordsB=ax_lower.transData,
        linestyle="--", color="grey", linewidth=0.8, zorder=1,
    )
    fig.add_artist(conn)


# -------------------------
# Plot: 2D
# -------------------------

def plot_whole_profile_2d(fig_x, fig_y,
                          cmap,
                          profiles: dict,
                          contexts: list,
                          outfile: str,
                          distance: int,
                          bins_st: int, bins_body: int, bins_en: int,
                          line_lw: float = 2.0,
                          ylim_dict: dict = None,
                          xpad: int = 10,
                          line_colors: list = None,
                          standard_distance_labels: bool = False):
    """
    2D version: nrows = len(contexts), one vertically stacked panel per context.
    """
    all_length, x_pos, x_labels, x_lines = make_xticks(
        distance, bins_st, bins_body, bins_en,
        standard_distance_labels=standard_distance_labels
    )

    nrows = len(contexts)
    labels = list(next(iter(profiles.values())).keys())  # sample labels used in the legend
    n_samples = len(labels)

    # Line colors: use user-defined colors first; otherwise derive colors from cmap
    if line_colors is not None:
        colors = line_colors
    else:
        raw_cmap = plt.get_cmap(cmap, max(3, n_samples))
        colors = [raw_cmap(i) for i in range(n_samples)]

    fig, axes = plt.subplots(nrows, 1, figsize=(fig_x, fig_y), sharex=True)
    if nrows == 1:
        axes = [axes]

    for i, cx in enumerate(contexts):
        ax = axes[i]
        ax.set_xlim(-xpad, all_length + xpad)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Curves
        for j, lab in enumerate(labels):
            y = profiles[cx][lab]["methylation"].values
            ax.plot(np.arange(len(y)), y, lw=line_lw, color=colors[j], label=lab, zorder=3)

        # TSS/TES vertical guide lines
        for xl in x_lines:
            ax.axvline(xl, linestyle="--", color="grey", linewidth=1, zorder=1)

        ax.set_ylabel(cx)

        if ylim_dict is not None:
            ax.set_ylim(*ylim_dict[cx])

        # Put the legend in the first panel
        if i == 0:
            ax.legend(frameon=False, loc="upper right")

    axes[-1].set_xticks(x_pos, x_labels)
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


# -------------------------
# Plot: 3D (layered)
# -------------------------

def plot_whole_profile_3d(fig_x, fig_y,
                          cmap,
                          profiles: dict,
                          contexts: list,
                          outfile: str,
                          distance: int,
                          bins_st: int, bins_body: int, bins_en: int,
                          edge: float = 0.10,
                          x_offset_pct: float = 20.0,
                          y_offset_pct: float = 35.0,
                          line_lw: float = 2.0,
                          ylim_dict: dict = None,
                          xpad: int = 10,
                          line_colors: list = None,
                          standard_distance_labels: bool = False):
    """
    3D version (pseudo-3D): each layer is a context (CG/CHG/CHH), and layering is achieved by offsetting the axes positions.

    - edge: figure margin (0-1, figure fraction)
    - x_offset_pct / y_offset_pct: offset percentage relative to the layer width/height (0-100)
      For example, x_offset_pct=20 means each lower layer shifts right by 0.2 * width
    """
    all_length, x_pos, x_labels, x_lines = make_xticks(
        distance, bins_st, bins_body, bins_en,
        standard_distance_labels=standard_distance_labels
    )

    n_layers = len(contexts)
    labels = list(next(iter(profiles.values())).keys())
    n_samples = len(labels)

    # Line colors: use user-defined colors first; otherwise derive colors from cmap
    if line_colors is not None:
        colors = line_colors
    else:
        raw_cmap = plt.get_cmap(cmap, max(3, n_samples))
        colors = [raw_cmap(i) for i in range(n_samples)]

    # Layer width and height: shrink each layer according to the offsets so that all layers fit within the canvas
    base_w = 1.0 - 2 * edge
    base_h = 1.0 - 2 * edge
    width = base_w / (1.0 + (x_offset_pct / 100.0) * max(0, n_layers - 1))
    height = base_h / (1.0 + (y_offset_pct / 100.0) * max(0, n_layers - 1))

    fig = plt.figure(figsize=(fig_x, fig_y))
    ax_list = []

    # Larger i means the layer is visually more in front (further right and lower)
    for i, cx in enumerate(contexts):
        left = edge + (x_offset_pct / 100.0) * i * width
        bottom = edge + (y_offset_pct / 100.0) * (n_layers - 1 - i) * height

        ax = fig.add_axes([left, bottom, width, height])
        ax_list.append(ax)

        # Transparent background so lower layers remain visible
        ax.set_facecolor("none")
        ax.patch.set_visible(False)

        ax.set_xlim(-xpad, all_length + xpad)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Show x ticks only on the front-most (lowest) layer
        if i != n_layers - 1:
            plt.setp(ax.get_xticklabels(), visible=False)

        # Curves
        for j, lab in enumerate(labels):
            y = profiles[cx][lab]["methylation"].values
            ax.plot(np.arange(len(y)), y, lw=line_lw, color=colors[j], label=lab, zorder=3)

        ax.set_ylabel(cx)
        ax.set_xticks(x_pos, x_labels)

        # TSS / TES vertical guide lines
        for xl in x_lines:
            ax.axvline(xl, linestyle="--", color="grey", linewidth=1, zorder=1)

        if ylim_dict is not None:
            ax.set_ylim(*ylim_dict[cx])

        # Legend: place it to the right of the front-most layer
        if i == n_layers - 1:
            ax.legend(frameon=False, loc=[1.02, 0.3])

    # zorder: ensure front layers are drawn above back layers
    for z, ax in enumerate(ax_list):
        ax.set_zorder(z)

    # Vertical dashed connectors between layers (TSS/TES)
    for i in range(len(ax_list) - 1):
        for xpos in x_lines:
            connect_vertical(xpos, ax_list[i], ax_list[i + 1])

    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


# -------------------------
# Main (CLI)
# -------------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Plot whole-profile methylation curves in 2D and/or 3D layered style.\n\n"
            "Input files are resolved by: <meth_dir>/<sample>_<context><whole_suffix>.\n"
            "Example: CX_gene/SRR9321764_CHH_profile.tsv (use --whole-suffix _profile.tsv)."
        )
    )

    # I/O
    ap.add_argument(
        "--meth-dir", default="gene",
        help=(
            "Directory containing methylation profile files (default: gene/). "
            "Files are resolved as: <meth_dir>/<sample>_<context><whole_suffix>. "
            "You can use CX_gene/ or gene/ etc."
        )
    )
    ap.add_argument(
        "--samples", required=True,
        help=(
            "Sample IDs used to build filenames, comma-separated. "
            "Example: SRR9321764,SRR8742373"
        )
    )
    ap.add_argument(
        "--labels", required=True,
        help=(
            "Display labels for legend, comma-separated and one-to-one with --samples. "
            "Example: WT,Mut"
        )
    )
    ap.add_argument(
        "--contexts", default="CG,CHG,CHH",
        help=(
            "Methylation contexts (comma-separated). This list also controls panel order in the output. "
            "Example: CG,CHG,CHH or CHH"
        )
    )
    ap.add_argument(
        "--whole-suffix", default="_whole_line.txt",
        help=(
            "Suffix of the whole-profile file. Default: _whole_line.txt. "
            "If your file is named like '<sample>_<context>_profile.tsv', set --whole-suffix _profile.tsv"
        )
    )
    ap.add_argument(
        "--out-prefix", default="whole_profile_methylation",
        help=(
            "Output prefix. The script writes '<prefix>.2D.pdf' and/or '<prefix>.3D.pdf' in the current folder."
        )
    )

    # Gene model / bins (used ONLY for ticks & vertical lines)
    ap.add_argument(
        "--distance", type=int, default=2000,
        help=(
            "Upstream/downstream window size in bp for axis labels (does NOT rescale your input). "
            "Used to annotate -distance/+distance on x-axis."
        )
    )
    ap.add_argument(
        "--bins-st", type=int, default=50,
        help=(
            "Number of bins for upstream segment (from -distance to TSS). "
            "Must match the binning used to generate your profile file."
        )
    )
    ap.add_argument(
        "--bins-body", type=int, default=100,
        help=(
            "Number of bins for gene body (TSS to TES). Must match your profile file."
        )
    )
    ap.add_argument(
        "--bins-en", type=int, default=50,
        help=(
            "Number of bins for downstream segment (TES to +distance). Must match your profile file."
        )
    )
    ap.add_argument(
        "--xpad", type=int, default=10,
        help=(
            "Extra padding (in bins) added to both ends of x-axis for aesthetics only."
        )
    )

    # Style
    ap.add_argument(
        "--cmap", default="tab10",
        help=(
            "Matplotlib colormap name used to assign colors when --line-colors is not provided. "
            "Examples: tab10, Set2, viridis, RdBu_r"
        )
    )
    ap.add_argument("--line-lw", type=float, default=2.0, help="Line width for curves.")
    ap.add_argument(
        "--line-colors", default=None,
        help=(
            "Optional custom line colors for samples, in the same order as --labels. "
            "Comma-separated; supports hex or named colors. "
            "Examples: '#FF3030,#00EE76' or 'red,blue'. "
            "If only 1 color is provided, all samples use it; if fewer than samples, colors will cycle."
        )
    )
    ap.add_argument("--fig-x", type=float, default=6.0, help="Figure width in inches.")
    ap.add_argument("--fig-y", type=float, default=8.0, help="Figure height in inches.")

    # 3D layout
    ap.add_argument(
        "--edge", type=float, default=0.10,
        help=(
            "Figure margin used by the 3D layered layout (0-1). Example: 0.10 means ~10%% blank around panels."
        )
    )
    ap.add_argument(
        "--x-offset", type=float, default=20.0,
        help=(
            "Horizontal offset per layer in 3D plot, as percent of layer width. "
            "Positive moves layers to the right."
        )
    )
    ap.add_argument(
        "--y-offset", type=float, default=35.0,
        help=(
            "Vertical offset per layer in 3D plot, as percent of layer height. "
            "Positive moves layers downward."
        )
    )

    # ylim control
    ap.add_argument(
        "--ylim", default=None,
        help=(
            "Optional y-axis limits. Format A (global): 'ymin,ymax'. "
            "Format B (per context in --contexts order): 'ymin,ymax;ymin,ymax;...'."
        )
    )

    # distance label style
    ap.add_argument(
        "--standard-distance-labels", action="store_true",
        help=(
            "Use standard labels: left '-distance' and right '+distance'. "
            "By default, the script keeps the legacy label style: left '+distance' and right '-distance'."
        )
    )

    # mode
    ap.add_argument(
        "--mode", default="both", choices=["2d", "3d", "both"],
        help="Which outputs to generate: 2d, 3d, or both (default: both)."
    )

    args = ap.parse_args()

    meth_dir = args.meth_dir
    samples = _split_csv(args.samples)
    labels = _split_csv(args.labels)
    contexts = _split_csv(args.contexts)

    ylim_dict = parse_ylim(args.ylim, contexts) if args.ylim else None
    user_line_colors = parse_line_colors(args.line_colors, len(labels))

    profiles = load_profiles(
        meth_dir=meth_dir,
        samples=samples,
        labels=labels,
        contexts=contexts,
        whole_suffix=args.whole_suffix
    )

    if args.mode in ("2d", "both"):
        plot_whole_profile_2d(
            fig_x=args.fig_x, fig_y=args.fig_y,
            cmap=args.cmap,
            profiles=profiles,
            contexts=contexts,
            outfile=f"{args.out_prefix}.2D.pdf",
            distance=args.distance,
            bins_st=args.bins_st, bins_body=args.bins_body, bins_en=args.bins_en,
            line_lw=args.line_lw,
            ylim_dict=ylim_dict,
            xpad=args.xpad,
            line_colors=user_line_colors,
            standard_distance_labels=args.standard_distance_labels
        )
    print(f"[OK] Saved: {args.out_prefix}.2D.pdf")

    if args.mode in ("3d", "both"):
        plot_whole_profile_3d(
            fig_x=args.fig_x, fig_y=args.fig_y,
            cmap=args.cmap,
            profiles=profiles,
            contexts=contexts,
            outfile=f"{args.out_prefix}.3D.pdf",
            distance=args.distance,
            bins_st=args.bins_st, bins_body=args.bins_body, bins_en=args.bins_en,
            edge=args.edge,
            x_offset_pct=args.x_offset,
            y_offset_pct=args.y_offset,
            line_lw=args.line_lw,
            ylim_dict=ylim_dict,
            xpad=args.xpad,
            line_colors=user_line_colors,
            standard_distance_labels=args.standard_distance_labels
        )
    print(f"[OK] Saved: {args.out_prefix}.3D.pdf")

if __name__ == "__main__":
    main()
