#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OmicsCanvas: Whole-profile plotter (2D + 3D stacked tracks)

This script plots 1D meta-profiles (mean/median across genes) from matrix files.

Input matrices are produced by:
  - omicscanvas_bam_to_gene_matrices.py  (recommended, standard naming)
  - 2caculate_gene_matrix_new_prefix.py (legacy)

Default output suffixes (standard naming, --naming standard):
  TSS  : *_tss_matrix.tsv
  gene : *_gene_profile_matrix.tsv
  TES  : *_tes_matrix.tsv

Legacy suffixes (--naming legacy):
  TSS  : *_start_matrix.txt
  gene : *_gene_body_matrix.txt
  TES  : *_end_matrix.txt

GROUP SPEC syntax (most important):
  - Comma (,) separates samples inside one panel (drawn as multiple lines).
  - Semicolon (;) separates panels within one COLUMN (stacked).
  - Pipe (|) starts a NEW COLUMN.

Example (2 columns; each column has 2 panels; each panel has 3 samples):
  --group   "G_Input,B_Input,R_Input;G_H3K4me3,B_H3K4me3,R_H3K4me3|G_H3K27me3,B_H3K27me3,R_H3K27me3;G_Pol2,B_Pol2,R_Pol2"   --names   "G,B,R;G,B,R|G,B,R;G,B,R"   --ylabels "Input;H3K4me3|H3K27me3;PolII"

File mapping rule (IMPORTANT):
  matrix_path = <matrix_dir>/<sample_prefix><suffix_by_gene_type>

Notes:
  - connect_vertical() is used as-is (no endpoint changes).
  - Run "python <script>.py -h" to see the full parameter list and defaults.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch


# ----------------------------
# Config loader (optional)
# ----------------------------
def load_config_defaults(config_path: str | None) -> dict:
    """
    Load a JSON config file and convert keys to argparse dest names.

    - Keys may use '-' or '_' (both accepted).
    - Some user-friendly aliases are supported (e.g., matrix_dir -> dir_chip).
    - Unknown keys are ignored with a warning (after the full parser is built).
    """
    if not config_path:
        return {}
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    cfg = json.loads(p.read_text(encoding='utf-8'))
    if not isinstance(cfg, dict):
        raise ValueError("--config must be a JSON object (dict).")

    norm = {}
    for k, v in cfg.items():
        kk = str(k).strip().lstrip('-').replace('-', '_')
        norm[kk] = v

    # Friendly aliases -> internal dest names
    alias = {
        'matrix_dir': 'dir_chip',
        'matrixdir': 'dir_chip',
        'input_dir': 'dir_chip',
        'inputdir': 'dir_chip',
        'out_file': 'out',
        'outfile': 'out',
        'figure_x': 'fig_x',
        'figure_y': 'fig_y',
    }
    out = {}
    for kk, v in norm.items():
        out[alias.get(kk, kk)] = v
    return out




# ----------------------------
# Parsing helpers
# ----------------------------
def _split_keep_empty(s: str, sep: str) -> list[str]:
    # consistent split with trimming (keeps empties removed)
    parts = [p.strip() for p in s.split(sep)]
    return [p for p in parts if p != ""]


def parse_columns(spec: str) -> list[list[list[str]]]:
    """
    Parse group-like spec into columns -> panels -> samples.
      columns separated by '|'
      panels separated by ';'
      samples separated by ','
    """
    cols: list[list[list[str]]] = []
    for col in _split_keep_empty(spec, "|"):
        panels: list[list[str]] = []
        for panel in _split_keep_empty(col, ";"):
            samples = _split_keep_empty(panel, ",")
            if samples:
                panels.append(samples)
        if panels:
            cols.append(panels)
    if not cols:
        raise ValueError("Spec parsed to empty. Check separators: ',' ';' '|'.")
    return cols


def flatten_panels(cols: list[list[list[str]]]) -> list[list[str]]:
    out: list[list[str]] = []
    for c in cols:
        out.extend(c)
    return out


def parse_panel_strings(spec: str, template_cols: list[list[list[str]]], what: str) -> list[list[list[str]]]:
    """
    Parse names/ylabels spec to match template columns/panels.
    For ylabels: each panel string may be single token (not comma split).
    For names  : each panel string is comma-split into sample labels.
    """
    # if user didn't put '|', treat as single stream over all panels
    if "|" not in spec and len(template_cols) > 1:
        # distribute sequentially across all panels
        all_panels = _split_keep_empty(spec, ";")
        flat = flatten_panels(template_cols)
        if len(all_panels) != len(flat):
            raise ValueError(f"{what} panels mismatch: need {len(flat)} panels but got {len(all_panels)}")
        # rebuild into columns
        idx = 0
        out_cols: list[list[list[str]]] = []
        for col in template_cols:
            out_col: list[list[str]] = []
            for _panel in col:
                s = all_panels[idx]
                idx += 1
                if what == "names":
                    out_col.append(_split_keep_empty(s, ","))
                else:  # ylabels
                    out_col.append([s])
            out_cols.append(out_col)
        return out_cols

    # user provided '|', parse per column
    cols_spec = _split_keep_empty(spec, "|")
    if len(cols_spec) != len(template_cols):
        raise ValueError(f"{what} column mismatch: need {len(template_cols)} columns but got {len(cols_spec)}")

    out_cols = []
    for col_str, col_tpl in zip(cols_spec, template_cols):
        panel_strs = _split_keep_empty(col_str, ";")
        if len(panel_strs) != len(col_tpl):
            raise ValueError(f"{what} panel mismatch in one column: need {len(col_tpl)} panels but got {len(panel_strs)}")
        out_col = []
        for s in panel_strs:
            if what == "names":
                out_col.append(_split_keep_empty(s, ","))
            else:
                out_col.append([s])
        out_cols.append(out_col)
    return out_cols


def enforce_names_match(group_cols: list[list[list[str]]], name_cols: list[list[list[str]]]) -> None:
    for ci, (gcol, ncol) in enumerate(zip(group_cols, name_cols)):
        for pi, (gpanel, npanel) in enumerate(zip(gcol, ncol)):
            if len(npanel) == 0:
                raise ValueError(f"names panel is empty at col={ci}, panel={pi}")
            if len(npanel) != len(gpanel):
                raise ValueError(
                    f"names count mismatch at col={ci}, panel={pi}: "
                    f"{len(gpanel)} samples in group, but {len(npanel)} in names"
                )


# ----------------------------
# Matrix I/O and profiling
# ----------------------------
def read_matrix(path: str, index_filter: str | None, scale: float) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    if index_filter:
        # default: keep transcript ".1"
        try:
            m = df.index.to_series().astype(str).str.contains(index_filter, regex=True)
            df = df.loc[m.values]
        except re.error as e:
            raise ValueError(f"Bad --index-filter regex: {index_filter}. {e}")
    if scale != 1.0:
        df = df * scale
    return df


def load_profiles_for_gene_type(
    dir_chip: str,
    group_cols: list[list[list[str]]],
    suffix: str,
    index_filter: str | None,
    scale: float,
) -> list[list[list[pd.DataFrame]]]:
    """
    Return same structure as group_cols, but each sample is a DataFrame (genes x bins).
    """
    out: list[list[list[pd.DataFrame]]] = []
    for col in group_cols:
        out_col: list[list[pd.DataFrame]] = []
        for panel in col:
            out_panel: list[pd.DataFrame] = []
            for sample_prefix in panel:
                fp = os.path.join(dir_chip, sample_prefix + suffix)
                if not os.path.exists(fp):
                    raise FileNotFoundError(f"Matrix not found: {fp}")
                out_panel.append(read_matrix(fp, index_filter=index_filter, scale=scale))
            out_col.append(out_panel)
        out.append(out_col)
    return out


def panel_profile(panel_dfs: list[pd.DataFrame], stat: str = "mean") -> list[np.ndarray]:
    """
    For each sample DataFrame (genes x bins), return 1D profile (bins,)
    """
    profiles: list[np.ndarray] = []
    for df in panel_dfs:
        if stat == "median":
            s = df.median(axis=0)
        else:
            s = df.mean(axis=0)
        profiles.append(s.to_numpy(dtype=float))
    return profiles


def get_y_limits_for_panel(panel_profiles: list[np.ndarray], ymin_zero: bool, pad_frac: float = 0.05) -> tuple[float, float]:
    arr = np.concatenate([p[np.isfinite(p)] for p in panel_profiles if p is not None and len(p) > 0])
    if arr.size == 0:
        return (0.0, 1.0)
    ymin = float(np.nanmin(arr))
    ymax = float(np.nanmax(arr))
    if ymin_zero:
        ymin = 0.0
    if ymax == ymin:
        ymax = ymin + 1.0
    pad = (ymax - ymin) * pad_frac
    return (ymin - pad, ymax + pad)


# ----------------------------
# X ticks / lines (inheritable from 2caculate params)
# ----------------------------
def build_x_guides(distance: int, bins_start: int, bins_gene_st: int, bins_gene_body: int, bins_gene_en: int, bins_end: int):
    xticks = {}
    xlines = {}

    # NOTE: Keep user's original label direction (they used "+distance" on left).
    xticks["TSS"] = {
        "pos": [0, bins_start / 2 - 1, bins_start - 1],
        "labels": ["+" + str(distance), "TSS", "-" + str(distance)],
    }
    xlines["TSS"] = [bins_start / 2 - 1]

    xticks["gene"] = {
        "pos": [0, bins_gene_st, bins_gene_st + bins_gene_body, bins_gene_st + bins_gene_body + bins_gene_en],
        "labels": ["+" + str(distance), "TSS", "TES", "-" + str(distance)],
    }
    xlines["gene"] = [bins_gene_st, bins_gene_st + bins_gene_body]

    xticks["TES"] = {
        "pos": [0, bins_end / 2 - 1, bins_end - 1],
        "labels": ["+" + str(distance), "TES", "-" + str(distance)],
    }
    xlines["TES"] = [bins_end / 2 - 1]

    return xticks, xlines


# ----------------------------
# Plot helpers
# ----------------------------
def connect_vertical(xpos, ax_upper, ax_lower):
    """
    Draw a vertical dashed line connecting 2 axes (coordinate transform).
    Do NOT change endpoint logic here (kept as user's convention).
    """
    y_bottom_of_upper = ax_upper.get_ylim()[0]
    y_top_of_lower = ax_lower.get_ylim()[0]

    fig = ax_upper.figure
    conn = ConnectionPatch(
        xyA=(xpos, y_bottom_of_upper), coordsA=ax_upper.transData,
        xyB=(xpos, y_top_of_lower), coordsB=ax_lower.transData,
        linestyle="--", color="grey", linewidth=0.8, zorder=1
    )
    conn.set_clip_on(False)
    fig.add_artist(conn)


def add_inside_xlines(ax, xs, *, color="grey", lw=0.8, ls="--", zorder=1):
    # x: data coords, y: axes coords (0..1)
    for x in xs:
        ax.plot([x, x], [0, 1], transform=ax.get_xaxis_transform(),
                color=color, lw=lw, ls=ls, zorder=zorder)


def _maybe_truncate_label(s: str, maxlen: int) -> str:
    if maxlen <= 0:
        return s
    if len(s) <= maxlen:
        return s
    # keep head+tail
    if maxlen <= 3:
        return s[:maxlen]
    return s[: maxlen - 1] + "…"




def calc_3d_layout(
    columns_panels,
    *,
    base_left=0.06,
    base_top=0.78,
    col_gap=0.06,
    col_width=0.38,
    panel_height=0.14,
    x_off=0.10,
    y_off=0.45,
    panel_gap=0.05,
    right_margin=0.02,
    bottom_margin=0.05,
    auto_scale=True,
):
    """
    Calculate nested rects[col_idx][panel_idx] = [left,bottom,width,height] for 3D stacked panels.

    columns_panels: list of columns; each column is a list of panels (panel order = top->bottom stacking).

    Key idea:
      - Each column's "effective width" is larger than col_width because deeper panels shift right by x_off.
        effective_width(col) = col_width * (1 + (n_panels-1)*x_off)
      - Column start x positions are computed cumulatively using effective widths, so columns never overlap.
      - If layout would overflow, auto-shrink col_width and/or panel_height (unless auto_scale=False).
    """
    if not columns_panels:
        raise ValueError("Empty columns_panels (nothing to plot).")

    ncols = len(columns_panels)
    n_panels_each = [len(col) for col in columns_panels]
    if any(n <= 0 for n in n_panels_each):
        raise ValueError("Each column must contain >=1 panel.")

    eps = 1e-9

    # -------------------------
    # Horizontal auto-scale
    # -------------------------
    # Total effective width = sum( (1+(n-1)*x_off) * col_width ) + (ncols-1)*col_gap
    eff_units = [1.0 + (n - 1) * float(x_off) for n in n_panels_each]
    avail_w = 1.0 - float(base_left) - float(right_margin) - float(col_gap) * (ncols - 1)
    if avail_w <= 0:
        raise ValueError("Horizontal space is negative. Reduce margins/gaps.")

    w_max = avail_w / sum(eff_units)
    used_w = float(col_width)
    if used_w > w_max:
        if not auto_scale:
            raise ValueError(
                f"3D layout overflow in X: need --col-width <= {w_max:.4f} (given {used_w:.4f}). "
                "Try smaller --col-width/--x-off, smaller --col-gap, or larger margins."
            )
        used_w = w_max

    # Column start positions using effective widths (prevents overlap/overflow)
    col_lefts = []
    cursor = float(base_left)
    for u in eff_units:
        col_lefts.append(cursor)
        cursor += u * used_w + float(col_gap)

    # -------------------------
    # Vertical auto-scale
    # -------------------------
    max_n = max(n_panels_each)

    # Total effective height for a column:
    # used_h + (n-1)*(y_off*used_h + panel_gap) = used_h*(1+(n-1)*y_off) + (n-1)*panel_gap
    avail_h = float(base_top) - float(bottom_margin)
    if avail_h <= 0:
        raise ValueError("Vertical space is negative. Increase --base-top or reduce --bottom-margin.")

    denom = 1.0 + (max_n - 1) * float(y_off)
    numer = avail_h - (max_n - 1) * float(panel_gap)
    h_max = numer / denom if denom > 0 else 0.0
    if h_max <= 0:
        raise ValueError(
            "Vertical space is insufficient for the requested stack.\n"
            "Try smaller --panel-gap/--y-off, smaller --panel-height, or increase --base-top."
        )

    used_h = float(panel_height)
    if used_h > h_max:
        if not auto_scale:
            raise ValueError(
                f"3D layout overflow in Y: need --panel-height <= {h_max:.4f} (given {used_h:.4f}). "
                "Try smaller --panel-height/--y-off/--panel-gap, or adjust margins."
            )
        used_h = h_max

    # -------------------------
    # Generate nested rects
    # -------------------------
    rects = []
    for col_idx, panels in enumerate(columns_panels):
        col_left = col_lefts[col_idx]
        col_rects = []
        for j, _panel in enumerate(panels):
            left = col_left + j * float(x_off) * used_w
            bottom = float(base_top) - used_h - j * (float(y_off) * used_h + float(panel_gap))
            rect = [left, bottom, used_w, used_h]

            # strict bounds check (should pass after auto-scale; allow epsilon)
            if (left < -eps or bottom < -eps or
                left + used_w > 1.0 + eps or bottom + used_h > 1.0 + eps):
                if not auto_scale:
                    raise ValueError(
                        f"3D layout overflow (col={col_idx}, panel={j}): {rect}\n"
                        "Try: smaller --col-width/--panel-height or smaller --x-off/--y-off/--panel-gap or larger margins."
                    )
                # Very conservative fallback shrink (one-shot)
                # shrink width to fit right boundary
                max_right = 1.0 - float(right_margin)
                max_w = (max_right - col_left) / (1.0 + j * float(x_off) + 1e-12)
                used_w2 = min(used_w, max_w)

                # shrink height to fit bottom boundary
                max_h = (float(base_top) - float(bottom_margin) - j * float(panel_gap)) / (1.0 + j * float(y_off) + 1e-12)
                used_h2 = min(used_h, max_h)

                used_w, used_h = max(used_w2, 1e-6), max(used_h2, 1e-6)
                left = col_left + j * float(x_off) * used_w
                bottom = float(base_top) - used_h - j * (float(y_off) * used_h + float(panel_gap))
                rect = [left, bottom, used_w, used_h]

            col_rects.append(rect)
        rects.append(col_rects)

    return rects, used_w, used_h

def plot_panel_axes_3d(fig, rect, profiles, colors, labels, ylabel, xguides, xlines, show_xticks, line_lw, ylim, ylabel_maxlen):
    ax = fig.add_axes(rect)
    ax.set_facecolor("none")
    ax.patch.set_visible(False)

    ax.set_xlim(-10, len(profiles[0]) + 10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for i, prof in enumerate(profiles):
        ax.plot(np.arange(len(prof)), prof, color=colors[i], lw=line_lw, label=labels[i], zorder=3)

    ax.set_ylabel(_maybe_truncate_label(ylabel, ylabel_maxlen))

    # vertical guides inside axis
    add_inside_xlines(ax, xlines, zorder=1)

    if show_xticks:
        ax.set_xticks(xguides["pos"])
        ax.set_xticklabels(xguides["labels"])
    else:
        ax.set_xticks(xguides["pos"])
        ax.set_xticklabels([])

    if ylim is not None:
        ax.set_ylim(*ylim)

    return ax


def plot_panel_axes_2d(fig, rect, profiles, colors, labels, ylabel, xguides, xlines, show_xticks, line_lw, ylim, ylabel_maxlen):
    ax = fig.add_axes(rect)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for i, prof in enumerate(profiles):
        ax.plot(np.arange(len(prof)), prof, color=colors[i], lw=line_lw, label=labels[i], zorder=3)

    ax.set_ylabel(_maybe_truncate_label(ylabel, ylabel_maxlen))

    for xl in xlines:
        ax.axvline(xl, ls="--", color="grey", lw=0.8, zorder=1)

    if show_xticks:
        ax.set_xticks(xguides["pos"])
        ax.set_xticklabels(xguides["labels"])
    else:
        ax.set_xticks(xguides["pos"])
        ax.set_xticklabels([])

    if ylim is not None:
        ax.set_ylim(*ylim)

    return ax


def build_colors(n: int, cmap: str | None, line_colors: str | None):
    if line_colors:
        cols = [c.strip() for c in line_colors.split(",") if c.strip()]
        if len(cols) < n:
            # cycle
            cols = (cols * (n // len(cols) + 1))[:n]
        return cols[:n]
    if cmap:
        cm = plt.get_cmap(cmap, n)
        return [cm(i) for i in range(n)]
    # fallback
    cm = plt.get_cmap("Set2", n)
    return [cm(i) for i in range(n)]


# ----------------------------
# Main plotting
# ----------------------------
def main():
    # ---- optional config (JSON) ----
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument(
        "--config",
        default=None,
        help=(
            "Path to a JSON config file. Keys may use '-' or '_' and will be mapped to argument names. "
            "CLI arguments always override config values."
        ),
    )
    pre_args, _ = pre.parse_known_args()
    cfg = load_config_defaults(pre_args.config)

    p = argparse.ArgumentParser(
        parents=[pre],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="OmicsCanvas: plot whole ChIP/RNA profiles from matrix files (2D or stacked 3D)."
    )
    p.add_argument(
        "--mode",
        choices=["2d", "3d"],
        required=True,
        help=(
            "Plot mode. '2d' = one axis per panel stacked vertically. "
            "'3d' = fake-3D stacking inside each column using axis offsets + vertical connectors."
        ),
    )

    p.add_argument(
        "--matrix-dir",
        dest="dir_chip",
        default="caculate_matrix",
        help=(
            "Directory that contains matrix files generated by omicscanvas_bam_to_gene_matrices.py (or the legacy generator). "
            "For each sample_prefix in --group, the script loads: <matrix-dir>/<sample_prefix><suffix>."
        ),
    )
    # Backward-compatible alias (hidden)
    p.add_argument(
        "--dir-chip",
        dest="dir_chip",
        help=argparse.SUPPRESS,
    )
    p.add_argument(
        "--group",
        required=True,
        help=(
            "Group specification string using separators: ',' samples in one panel; ';' panels in a column; '|' new column. "
            "Each token is a sample_prefix used to build file path: <matrix-dir>/<sample_prefix><suffix>."
        ),
    )
    p.add_argument(
        "--names",
        required=True,
        help=(
            "Legend text for each sample line. Structure MUST match --group (same columns/panels/sample-count). "
            "Use the same separators ',', ';', '|'. Example: 'G,B,R;G,B,R|G,B,R;G,B,R'."
        ),
    )
    p.add_argument(
        "--ylabels",
        required=True,
        help=(
            "Y-axis label per panel. Panels must match --group panels. Use ';' and optional '|' (no comma split). "
            "Example: 'Input;H3K4me3|H3K27me3;PolII'."
        ),
    )
    p.add_argument(
        "--gene-type",
        choices=["TSS", "gene", "TES"],
        default="gene",
        help=(
            "Which region to plot. This determines the x-axis guide (TSS/TES positions) and which matrix file is loaded. "
            "Use --naming to select the standard vs legacy matrix filename suffixes."
        ),
    )

    p.add_argument(
        "--naming",
        choices=["standard", "legacy"],
        default="standard",
        help=(
            "Matrix filename suffix convention. "
            "standard expects *_tss_matrix.tsv / *_gene_profile_matrix.tsv / *_tes_matrix.tsv. "
            "legacy expects *_start_matrix.txt / *_gene_body_matrix.txt / *_end_matrix.txt. "
            "Advanced: you can override suffixes via hidden --suffix-tss/--suffix-gene/--suffix-tes options."
        ),
    )

    # Inherited parameters (MUST match the matrix generation script)
    p.add_argument(
        "--distance",
        type=int,
        default=2000,
        help=(
            "Flanking distance (bp) used ONLY for x tick labels. "
            "This should be the same 'distance' you used when generating matrices."
        ),
    )
    p.add_argument(
        "--bins-start",
        type=int,
        default=100,
        help=(
            "Total number of bins in the TSS window matrix. "
            "Used to place TSS tick and the vertical guide line for gene-type=TSS."
        ),
    )
    p.add_argument(
        "--bins-gene-st",
        type=int,
        default=50,
        help=(
            "Number of bins for the upstream/TSS part inside the gene-profile matrix. "
            "Used to place the 'TSS' guide line when gene-type=gene."
        ),
    )
    p.add_argument(
        "--bins-gene-body",
        type=int,
        default=100,
        help=(
            "Number of bins representing the scaled gene body inside the gene-profile matrix. "
            "Used to place the 'TES' guide line when gene-type=gene."
        ),
    )
    p.add_argument(
        "--bins-gene-en",
        type=int,
        default=50,
        help=(
            "Number of bins for the downstream/TES-flank part inside the gene-profile matrix. "
            "Total bins for gene-type=gene is bins-gene-st + bins-gene-body + bins-gene-en."
        ),
    )
    p.add_argument(
        "--bins-end",
        type=int,
        default=100,
        help=(
            "Total number of bins in the TES window matrix. "
            "Used to place TES tick and the vertical guide line for gene-type=TES."
        ),
    )

    # File suffix mapping (sample_prefix + suffix -> matrix filename)
    # Advanced/compatibility: you can still override suffixes explicitly (hidden in -h).
    p.add_argument("--suffix-tss", "--suffix-start", dest="suffix_start", default="_tss_matrix.tsv", help=argparse.SUPPRESS)
    p.add_argument("--suffix-gene", dest="suffix_gene", default="_gene_profile_matrix.tsv", help=argparse.SUPPRESS)
    p.add_argument("--suffix-tes", "--suffix-end", dest="suffix_end", default="_tes_matrix.tsv", help=argparse.SUPPRESS)

    # Data options
    p.add_argument(
        "--index-filter",
        default=r"\.1$",
        help=(
            "Regex filter applied to the FIRST COLUMN index (gene/transcript IDs) of the matrix. "
            "Default keeps only '.1' transcripts (common for isoform filtering). "
            "Use empty string '' to disable filtering."
        ),
    )
    p.add_argument(
        "--scale",
        type=float,
        default=10.0,
        help=(
            "Multiply matrix values by this factor AFTER reading. "
            "Use this to convert units (e.g., 0.1 -> 1) or to match your plotting convention."
        ),
    )
    p.add_argument(
        "--stat",
        choices=["mean", "median"],
        default="mean",
        help=(
            "How to aggregate across rows (genes) to produce a 1D profile per sample: mean or median."
        ),
    )

    # Style
    p.add_argument(
        "--fig-x",
        type=float,
        default=10.0,
        help="Figure width in inches (matplotlib figsize[0]).",
    )
    p.add_argument(
        "--fig-y",
        type=float,
        default=10.0,
        help="Figure height in inches (matplotlib figsize[1]).",
    )
    p.add_argument("--cmap", default="Set2", help="Colormap to generate line colors (ignored if --line-colors given).")
    p.add_argument("--line-colors", default=None, help="Comma colors for lines, e.g. 'forestgreen,dodgerblue,orangered'.")
    p.add_argument("--line-lw", type=float, default=1.0, help="Line width for each sample profile curve.")
    p.add_argument("--ylabel-maxlen", type=int, default=18, help="Truncate long ylabels. 0 to disable.")
    p.add_argument("--legend", action="store_true", help="Show legend on the first panel of each column.")
    p.add_argument("--legend-loc", default="outside", choices=["outside", "inside"], help="Legend placement.")
    p.add_argument(
        "--ylim",
        default=None,
        help=(
            "Manual y-axis limits applied to ALL panels, formatted as 'ymin,ymax'. "
            "If set, this overrides --ylims and auto-scaling."
        ),
    )
    p.add_argument(
        "--ylims",
        default=None,
        help=(
            "Manual y-axis limits PER panel. Use ';' (and optional '|') to match panel structure. "
            "Each panel token is 'ymin,ymax'. Example: '0,1;0,2|0,1;0,3'. "
            "If set, this overrides auto-scaling but does NOT override --ylim (global)."
        ),
    )
    p.add_argument("--share-y", default="none", choices=["none", "col", "all"], help="Share y-limits within a column or across all panels.")
    p.add_argument("--ymin-zero", action="store_true", help="Force panel ymin to 0 when auto-scaling y.")

    # 3D layout (ONLY used when --mode 3d)
    # These are all in *figure fraction* coordinates (0..1), not inches.
    p.add_argument(
        "--base-left",
        type=float,
        default=0.07,
        help=(
            "Left margin (figure fraction) where the FIRST column starts. "
            "Smaller -> more space for columns; larger -> more space for ylabels."
        ),
    )
    p.add_argument(
        "--base-top",
        type=float,
        default=0.93,
        help=(
            "Top anchor (figure fraction) for the FIRST (top) panel in each column. "
            "Decrease if titles/legend overflow; increase to give more vertical room."
        ),
    )
    p.add_argument(
        "--col-width",
        type=float,
        default=0.38,
        help=(
            "Width of the front-most panel in each column (figure fraction). "
            "If you have many columns or large x-off, reduce this."
        ),
    )
    p.add_argument(
        "--panel-height",
        type=float,
        default=0.14,
        help=(
            "Height of each panel (figure fraction). If you have many stacked panels, reduce this."
        ),
    )
    p.add_argument(
        "--x-off",
        type=float,
        default=0.10,
        help=(
            "Horizontal offset multiplier for deeper panels in the same column. "
            "Panel j is shifted by j * x_off * col_width. Larger -> stronger '3D' depth but consumes width."
        ),
    )
    p.add_argument(
        "--y-off",
        type=float,
        default=0.45,
        help=(
            "Vertical offset multiplier for deeper panels in the same column. "
            "Panel j is shifted downward by j * y_off * panel_height (plus panel-gap)."
        ),
    )
    p.add_argument(
        "--col-gap",
        type=float,
        default=0.06,
        help="Gap between columns (figure fraction). Increase to avoid overlap between adjacent columns.",
    )
    p.add_argument(
        "--panel-gap",
        type=float,
        default=0.05,
        help="Extra vertical gap between panels inside a column (figure fraction).",
    )

    p.add_argument(
        "--strict-layout",
        action="store_true",
        help=(
            "3D mode only. Do NOT auto-rescale the layout if panels/columns overflow the figure; instead raise an error. "
            "This is useful when you want reproducible geometry."
        ),
    )

    # Backward-compatible flags (hidden)
    p.add_argument("--auto-scale", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--no-auto-scale", action="store_true", help=argparse.SUPPRESS)

    # Output
    p.add_argument(
        "--out",
        required=True,
        help="Output figure path. Extension controls format (e.g., .pdf, .png, .svg).",
    )


    # Apply config defaults (ignore unknown keys)
    known = {a.dest for a in p._actions}
    unknown = sorted([k for k in cfg.keys() if k not in known])
    if unknown:
        print(f"[WARN] Unknown keys in --config ignored: {', '.join(unknown)}", file=sys.stderr)
    p.set_defaults(**{k: v for k, v in cfg.items() if k in known})

    args = p.parse_args()

    # ---- resolve suffix convention & layout strictness ----
    argv = set(sys.argv[1:])

    # Layout scaling: ON by default; can be disabled by --strict-layout (preferred) or legacy --no-auto-scale
    args.auto_scale = True
    if getattr(args, 'strict_layout', False) or getattr(args, 'no_auto_scale', False):
        args.auto_scale = False

    # Suffix convention (standard vs legacy), unless user explicitly overrides hidden suffix args
    user_set_start = ('--suffix-start' in argv) or ('--suffix-tss' in argv)
    user_set_gene  = '--suffix-gene'  in argv
    user_set_end   = ('--suffix-end' in argv) or ('--suffix-tes' in argv)

    if getattr(args, 'naming', 'standard') == 'legacy':
        if not user_set_start: args.suffix_start = '_start_matrix.txt'
        if not user_set_gene:  args.suffix_gene  = '_gene_body_matrix.txt'
        if not user_set_end:   args.suffix_end   = '_end_matrix.txt'
    else:
        if not user_set_start: args.suffix_start = '_tss_matrix.tsv'
        if not user_set_gene:  args.suffix_gene  = '_gene_profile_matrix.tsv'
        if not user_set_end:   args.suffix_end   = '_tes_matrix.tsv'

    dir_chip = args.dir_chip
    if not dir_chip.endswith("/"):
        dir_chip += "/"

    group_cols = parse_columns(args.group)
    name_cols = parse_panel_strings(args.names, group_cols, what="names")
    ylabel_cols = parse_panel_strings(args.ylabels, group_cols, what="ylabels")
    enforce_names_match(group_cols, name_cols)

    # suffix by gene type
    suffix = {"TSS": args.suffix_start, "gene": args.suffix_gene, "TES": args.suffix_end}[args.gene_type]

    index_filter = args.index_filter if args.index_filter.strip() != "" else None

    xticks_dict, xlines_dict = build_x_guides(
        distance=args.distance,
        bins_start=args.bins_start,
        bins_gene_st=args.bins_gene_st,
        bins_gene_body=args.bins_gene_body,
        bins_gene_en=args.bins_gene_en,
        bins_end=args.bins_end
    )
    xguides = xticks_dict[args.gene_type]
    xlines = xlines_dict[args.gene_type]

    # Load matrices for the selected gene_type
    data_cols = load_profiles_for_gene_type(
        dir_chip=dir_chip,
        group_cols=group_cols,
        suffix=suffix,
        index_filter=index_filter,
        scale=args.scale
    )

    # Build profiles + ylims
    # panel_profiles[col][panel] = list[np.ndarray] where list length = n_samples
    panel_profiles: list[list[list[np.ndarray]]] = []
    for col in data_cols:
        col_profiles = []
        for panel_dfs in col:
            col_profiles.append(panel_profile(panel_dfs, stat=args.stat))
        panel_profiles.append(col_profiles)

    # Prepare ylims:
    manual_all = None
    if args.ylim:
        a, b = [float(x.strip()) for x in args.ylim.split(",")]
        manual_all = (a, b)

    manual_per_panel = None
    if args.ylims:
        # parse like ylabels
        ycols_raw = parse_panel_strings(args.ylims, group_cols, what="ylabels")
        # each panel raw like ['0,1']
        manual_per_panel = []
        for col in ycols_raw:
            mcol = []
            for panel in col:
                s = panel[0]
                a, b = [float(x.strip()) for x in s.split(",")]
                mcol.append((a, b))
            manual_per_panel.append(mcol)

    # auto ylims per panel
    auto_ylims: list[list[tuple[float, float]]] = []
    for ci, col in enumerate(panel_profiles):
        col_ylims = []
        for pi, profs in enumerate(col):
            col_ylims.append(get_y_limits_for_panel(profs, ymin_zero=args.ymin_zero))
        auto_ylims.append(col_ylims)

    # apply sharing
    final_ylims = [[None for _ in col] for col in panel_profiles]

    if manual_all is not None:
        for ci in range(len(final_ylims)):
            for pi in range(len(final_ylims[ci])):
                final_ylims[ci][pi] = manual_all
    elif manual_per_panel is not None:
        for ci in range(len(final_ylims)):
            for pi in range(len(final_ylims[ci])):
                final_ylims[ci][pi] = manual_per_panel[ci][pi]
    else:
        if args.share_y == "all":
            ymin = min(auto_ylims[ci][pi][0] for ci in range(len(auto_ylims)) for pi in range(len(auto_ylims[ci])))
            ymax = max(auto_ylims[ci][pi][1] for ci in range(len(auto_ylims)) for pi in range(len(auto_ylims[ci])))
            for ci in range(len(final_ylims)):
                for pi in range(len(final_ylims[ci])):
                    final_ylims[ci][pi] = (ymin, ymax)
        elif args.share_y == "col":
            for ci in range(len(final_ylims)):
                ymin = min(y[0] for y in auto_ylims[ci])
                ymax = max(y[1] for y in auto_ylims[ci])
                for pi in range(len(final_ylims[ci])):
                    final_ylims[ci][pi] = (ymin, ymax)
        else:
            for ci in range(len(final_ylims)):
                for pi in range(len(final_ylims[ci])):
                    final_ylims[ci][pi] = auto_ylims[ci][pi]

    # Colors per panel (based on number of samples in that panel; usually constant within dataset)
    # We will build per panel to allow different sample counts across panels if needed.
    panel_colors: list[list[list]] = []
    for ci, col in enumerate(group_cols):
        col_colors = []
        for pi, panel in enumerate(col):
            col_colors.append(build_colors(len(panel), args.cmap, args.line_colors))
        panel_colors.append(col_colors)

    # ---------------- plot ----------------
    fig = plt.figure(figsize=(args.fig_x, args.fig_y))

    if args.mode == "2d":
        # Simple vertical stacking across ALL panels (flatten columns left->right)
        flat_panels = []
        flat_names = []
        flat_ylabels = []
        flat_ylims = []
        flat_colors = []

        for ci, col in enumerate(panel_profiles):
            for pi, profs in enumerate(col):
                flat_panels.append((ci, pi, profs))
                flat_names.append(name_cols[ci][pi])
                flat_ylabels.append(ylabel_cols[ci][pi][0])
                flat_ylims.append(final_ylims[ci][pi])
                flat_colors.append(panel_colors[ci][pi])

        n = len(flat_panels)
        top = 0.93
        bottom = 0.08
        left = 0.10
        width = 0.78
        gap = 0.02
        h = (top - bottom - gap * (n - 1)) / n
        if h <= 0:
            raise ValueError("Not enough vertical space for 2D panels. Increase fig_y or reduce panel count.")

        ax_list = []
        for i, (_ci, _pi, profs) in enumerate(flat_panels):
            rect = [left, top - (i + 1) * h - i * gap, width, h]
            show_xt = (i == n - 1)
            ax = plot_panel_axes_2d(
                fig, rect, profs, flat_colors[i], flat_names[i], flat_ylabels[i],
                xguides, xlines, show_xticks=show_xt, line_lw=args.line_lw,
                ylim=flat_ylims[i], ylabel_maxlen=args.ylabel_maxlen
            )
            ax_list.append(ax)

        if args.legend and ax_list:
            if args.legend_loc == "outside":
                ax_list[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
            else:
                ax_list[0].legend(loc="upper right", frameon=False)

    else:
        # 3D "stacked" within each column using offsets
        rects, used_w, used_h = calc_3d_layout(
            group_cols,
            base_left=args.base_left,
            base_top=args.base_top,
            col_width=args.col_width,
            panel_height=args.panel_height,
            x_off=args.x_off,
            y_off=args.y_off,
            col_gap=args.col_gap,
            panel_gap=args.panel_gap,
            auto_scale=args.auto_scale
        )

        # plot per column
        for ci, col in enumerate(panel_profiles):
            ax_list_col = []
            for pi, profs in enumerate(col):
                rect = rects[ci][pi]
                show_xt = (pi == len(col) - 1)
                ax = plot_panel_axes_3d(
                    fig, rect, profs,
                    panel_colors[ci][pi],
                    name_cols[ci][pi],
                    ylabel_cols[ci][pi][0],
                    xguides, xlines,
                    show_xticks=show_xt,
                    line_lw=args.line_lw,
                    ylim=final_ylims[ci][pi],
                    ylabel_maxlen=args.ylabel_maxlen
                )
                ax_list_col.append(ax)

            # zorder: lower panels under upper panels (like your manual reverse)
            for z, ax in enumerate(reversed(ax_list_col)):
                ax.set_zorder(z + 1)

            # connectors between stacked panels inside this column
            for j in range(len(ax_list_col) - 1):
                for xpos in xlines:
                    connect_vertical(xpos, ax_list_col[j], ax_list_col[j + 1])

            # legend once per column
            if args.legend and ax_list_col:
                if args.legend_loc == "outside":
                    ax_list_col[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
                else:
                    ax_list_col[0].legend(loc="upper right", frameon=False)

    # Save
    out = args.out
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Saved: {out}")


if __name__ == "__main__":
    main()
