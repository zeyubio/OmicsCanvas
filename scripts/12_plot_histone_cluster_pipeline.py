#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""omicscanvas_histone_cluster_pipeline.py

A practical analysis pipeline for ChIP-seq matrices produced by
`omicscanvas_bam_to_gene_matrices.py` (formerly `2caculate_gene_matrix_new_prefix.py`).

Key updates vs v2
-----------------
1) Expression (FPKM) plot not empty:
   - Adds robust ID normalization so matrix IDs (often transcript IDs like ...T00001.1)
     can be matched to expression IDs (often gene IDs like ...G00001).

2) Huge panel heatmaps:
   - Outputs PNG instead of PDF (PDF becomes enormous).
   - Produces TWO PNGs per region:
       (a) annotated: titles + xticks + boundary lines
       (b) clean: no titles/ticks/lines (easy for others to edit)

Inputs
------
1) ChIP-seq matrices (TSV): index column is ID, columns are bin_1..bin_N.
   Naming: <prefix> + <suffix>
     - TSS:  <prefix> + _tss_matrix.tsv
     - gene: <prefix> + _gene_profile_matrix.tsv
     - TES:  <prefix> + _tes_matrix.tsv

   Legacy suffixes are still supported via --naming legacy or --suffix-*.

2) FPKM table (TSV): index_col=0 is gene ID, replicates in columns.

Pipeline
--------
- Select marks and treatments using three compact strings:
    --in-group   : prefixes per mark, blocks separated by ';', within block by ','
    --in-ylabels : mark names, separated by ',' (must align with in_group blocks)
    --in-names   : treatment names per block, separated by ';' (align with in_group)

- Load matrices for selected regions (TSS/gene/TES).
- Normalize IDs to align matrices with expression (auto by default).
- For the CLUSTER region:
    per-gene zscore (axis=1) for each track (mark x treatment)
    -> concatenate tracks -> feature matrix -> KMeans clustering

Outputs
-------
1) Panel heatmaps (PNG) for each requested region:
   - *_panel_heatmaps_annot.png
   - *_panel_heatmaps_clean.png

2) Per-cluster mean profiles (PDF; compact; OK in PDF).
3) Expression boxplot per cluster (PDF).

Example
-------
python omicscanvas_histone_cluster_pipeline.py \
  --matrix-dir caculate_matrix \
  --in-group "G_H3K4me1,B_H3K4me1,R_H3K4me1;G_H3K4me3,B_H3K4me3,R_H3K4me3" \
  --in-ylabels "H3K4me1,H3K4me3" \
  --in-names "G,B,R;G,B,R" \
  --cluster-region gene \
  --plot-regions TSS,gene,TES \
  --fpkm FPKM_all.txt \
  --expr-cols "G=0,1,2;B=3,4,5;R=6,7,8" \
  --k 8 \
  --out-prefix G_B_R \
  --outdir out_histone

"""

# Headless-safe Matplotlib setup
# - Force a non-interactive backend (Agg)
# - Ensure MPLCONFIGDIR is writable (Matplotlib may crash/hang otherwise)

import os

# Put Matplotlib cache next to this script (typically under /mnt/data),
# not in $HOME (often read-only) and not in current working directory (may be '/').
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('MPLCONFIGDIR', os.path.join(_SCRIPT_DIR, '.mplconfig'))

import matplotlib
matplotlib.use('Agg')

import argparse
from pathlib import Path
import gzip
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Callable

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import zscore
from sklearn.cluster import KMeans
from matplotlib.backends.backend_pdf import PdfPages


# -----------------------------
# Parsing helpers
# -----------------------------

def _split_nonempty(s: str, sep: str) -> List[str]:
    return [x.strip() for x in s.split(sep) if x.strip()]


def parse_in_group(in_group: str, in_ylabels: str, in_names: str) -> Tuple[List[str], List[List[str]], List[List[str]]]:
    """Parse your compact strings.

    Returns
    -------
    marks: list[str]
    prefixes_per_mark: list[list[str]]
    treatments_per_mark: list[list[str]]
    """
    prefixes_groups = _split_nonempty(in_group, ';')
    marks = _split_nonempty(in_ylabels, ',')
    names_groups = _split_nonempty(in_names, ';')

    if len(prefixes_groups) != len(marks):
        raise ValueError(f"in_group blocks ({len(prefixes_groups)}) != in_ylabels marks ({len(marks)})")
    if len(names_groups) != len(prefixes_groups):
        raise ValueError(f"in_names blocks ({len(names_groups)}) != in_group blocks ({len(prefixes_groups)})")

    prefixes_per_mark: List[List[str]] = []
    treatments_per_mark: List[List[str]] = []

    for i in range(len(prefixes_groups)):
        prefixes = _split_nonempty(prefixes_groups[i], ',')
        treatments = _split_nonempty(names_groups[i], ',')
        if len(prefixes) != len(treatments):
            raise ValueError(f"mark[{i}] '{marks[i]}' prefixes ({len(prefixes)}) != treatments ({len(treatments)})")
        prefixes_per_mark.append(prefixes)
        treatments_per_mark.append(treatments)

    return marks, prefixes_per_mark, treatments_per_mark


def parse_expr_cols(expr_cols: str) -> Dict[str, List[int]]:
    """Parse expression column groups.

    Input format: "G=0,1,2;B=3,4,5" (0-based indices).
    """
    out: Dict[str, List[int]] = {}
    for part in _split_nonempty(expr_cols, ';'):
        if '=' not in part:
            raise ValueError(f"Bad --expr-cols item: {part}")
        k, v = part.split('=', 1)
        cols = [int(x) for x in _split_nonempty(v, ',')]
        out[k.strip()] = cols
    if not out:
        raise ValueError("--expr-cols parsed empty")
    return out


# -----------------------------
# ID normalization
# -----------------------------

def _id_identity(x: str) -> str:
    return x


def _id_drop_isoform(x: str) -> str:
    # remove trailing ".<digits>" (e.g., ".1")
    return re.sub(r"\.\d+$", "", x)


def _id_tx_to_gene(x: str) -> str:
    """Common plant annotation style:

    transcript:  Ca59v2Chr01T00001.1
    gene:        Ca59v2Chr01G00001

    Rule: replace the last 'T<digits>(.isoform)?' with 'G<digits>'.
    If pattern not matched, fall back to drop isoform only.
    """
    m = re.match(r"^(.*)T(\d+)(?:\.\d+)?$", x)
    if m:
        return f"{m.group(1)}G{m.group(2)}"
    return _id_drop_isoform(x)


def _make_id_transform(mode: str, regex_pat: str = "", regex_repl: str = "") -> Callable[[str], str]:
    if mode == 'none':
        return _id_identity
    if mode == 'drop_isoform':
        return _id_drop_isoform
    if mode == 'tx_to_gene':
        return _id_tx_to_gene
    if mode == 'regex':
        if not regex_pat:
            raise ValueError("--id-mode regex requires --id-regex")
        pat = re.compile(regex_pat)

        def _f(x: str) -> str:
            return pat.sub(regex_repl, x)

        return _f
    raise ValueError(f"Unknown id mode: {mode}")


def pick_best_id_mode(matrix_ids: List[str], expr_ids: List[str], allow: List[str]) -> str:
    """Pick a matrix->expr ID normalization mode that maximizes overlap."""
    expr_set = set(expr_ids)
    best_mode = allow[0]
    best_overlap = -1

    for mode in allow:
        f = _make_id_transform(mode)
        mapped = [f(x) for x in matrix_ids]
        overlap = len(set(mapped) & expr_set)
        if overlap > best_overlap:
            best_overlap = overlap
            best_mode = mode

    return best_mode


def apply_id_transform_and_collapse(df: pd.DataFrame,
                                   id_func: Callable[[str], str],
                                   collapse: str = 'mean') -> pd.DataFrame:
    """Transform index IDs -> new IDs and collapse duplicates."""
    new_ids = pd.Index([id_func(x) for x in df.index.astype(str)], name='ID')
    df2 = df.copy()
    df2.index = new_ids

    if df2.index.has_duplicates:
        if collapse == 'mean':
            df2 = df2.groupby(df2.index).mean(numeric_only=True)
        elif collapse == 'sum':
            df2 = df2.groupby(df2.index).sum(numeric_only=True)
        elif collapse == 'first':
            df2 = df2[~df2.index.duplicated(keep='first')]
        else:
            raise ValueError(f"Unknown --collapse-duplicates: {collapse}")

    df2.index = df2.index.astype(str)
    df2.index.name = 'ID'
    return df2


# -----------------------------
# I/O and transforms
# -----------------------------

def read_matrix(path: str,
                scale_factor: float = 10.0,
                isoform_suffix: str = ".1",
                bin_prefix: str = "bin_") -> pd.DataFrame:
    """Read a matrix TSV into gene x bins DataFrame."""
    df = pd.read_csv(path, sep='\t', index_col=0)
    df.index = df.index.astype(str)
    df.index.name = 'ID'

    # keep only one isoform if requested (common in your earlier scripts)
    if isoform_suffix:
        keep = [g for g in df.index if g.endswith(isoform_suffix)]
        df = df.loc[keep].copy()

    df = df.apply(pd.to_numeric, errors='coerce')

    bin_cols = [c for c in df.columns if str(c).startswith(bin_prefix)]
    if bin_cols:
        bin_cols = sorted(bin_cols, key=lambda x: int(str(x).split('_')[1]))
        df = df[bin_cols]

    if scale_factor is not None:
        df = df * float(scale_factor)

    return df


def zscore_by_gene(df: pd.DataFrame) -> pd.DataFrame:
    arr = zscore(df.values, axis=1, ddof=0, nan_policy='omit')
    return pd.DataFrame(arr, index=df.index, columns=df.columns)


def gene_intersection(dfs: List[pd.DataFrame]) -> List[str]:
    if not dfs:
        raise ValueError("No matrices provided to compute gene intersection")
    sets = [set(d.index) for d in dfs]
    common = set.intersection(*sets)
    return sorted(common)


# -----------------------------
# Plotting helpers
# -----------------------------

@dataclass
class RegionSpec:
    name: str
    bins_total: int
    xticks_pos: List[int]
    xticks_labels: List[str]
    vlines_pos: List[int]


def build_region_spec(region: str,
                      distance: int,
                      bins_start: int,
                      bins_end: int,
                      bins_gene_st: int,
                      bins_gene_body: int,
                      bins_gene_en: int) -> RegionSpec:
    """Define xticks and boundary lines."""
    if region == 'TSS':
        mid = bins_start // 2
        return RegionSpec(
            name='TSS',
            bins_total=bins_start,
            xticks_pos=[0, mid - 1, bins_start - 1],
            xticks_labels=[f"-{distance}", "TSS", f"+{distance}"],
            vlines_pos=[mid - 1],
        )
    if region == 'TES':
        mid = bins_end // 2
        return RegionSpec(
            name='TES',
            bins_total=bins_end,
            xticks_pos=[0, mid - 1, bins_end - 1],
            xticks_labels=[f"-{distance}", "TES", f"+{distance}"],
            vlines_pos=[mid - 1],
        )
    if region == 'gene':
        total = bins_gene_st + bins_gene_body + bins_gene_en
        return RegionSpec(
            name='gene',
            bins_total=total,
            xticks_pos=[0, bins_gene_st - 1, bins_gene_st + bins_gene_body - 1, total - 1],
            xticks_labels=[f"-{distance}", "TSS", "TES", f"+{distance}"],
            vlines_pos=[bins_gene_st, bins_gene_st + bins_gene_body],
        )
    raise ValueError(f"Unknown region: {region} (use TSS/gene/TES)")


def _parse_figsize_str(s: str):
    """Parse "W,H" (inches) -> (W, H) or None."""
    s = (s or '').strip()
    if not s:
        return None
    parts = [p.strip() for p in s.split(',') if p.strip()]
    if len(parts) != 2:
        raise ValueError(f"Bad --panel-figsize: {s!r}. Use 'W,H' e.g. '30,15'.")
    return float(parts[0]), float(parts[1])


def _sanitize_filename(s: str) -> str:
    # Keep letters, numbers, dot, dash, underscore; replace others with '_'
    return ''.join(ch if ch.isalnum() or ch in '._-+' else '_' for ch in s)


def plot_panel_heatmaps_png(track_mats_z: List[Tuple[str, pd.DataFrame]],
                            labels: pd.Series,
                            region_spec: RegionSpec,
                            out_png_annot: str,
                            out_png_clean: str,
                            cmap: str = 'RdBu_r',
                            vmin: float = -2,
                            vmax: float = 2,
                            panel_cols: Optional[int] = None,
                            panel_rows: Optional[int] = None,
                            cell_w: float = 3.5,
                            cell_h: float = 4.0,
                            figsize: Optional[Tuple[float, float]] = None,
                            dpi: int = 200):
    """Small-multiple heatmaps (one per track) -> TWO PNGs.

    - annotated: titles + xticks + boundary lines
    - clean    : no titles/ticks/lines (easy editing)

    Parameters
    ----------
    panel_cols / panel_rows
        Control grid layout. If only rows given, cols inferred.
    cell_w / cell_h
        Per-panel size (inches) used when figsize not provided.
    figsize
        Explicit overall figure size (inches). Overrides cell_w/cell_h.
    """

    order = labels.sort_values().index
    n_tracks = len(track_mats_z)
    if n_tracks == 0:
        raise ValueError('track_mats_z is empty')

    # Determine grid
    ncols = panel_cols
    nrows = panel_rows
    if ncols is None and nrows is None:
        ncols = min(9, n_tracks)
        nrows = int(np.ceil(n_tracks / ncols))
    elif ncols is None and nrows is not None:
        ncols = int(np.ceil(n_tracks / nrows))
    elif ncols is not None and nrows is None:
        nrows = int(np.ceil(n_tracks / ncols))

    if ncols <= 0 or nrows <= 0:
        raise ValueError(f'Invalid panel grid: rows={nrows}, cols={ncols}')
    if nrows * ncols < n_tracks:
        raise ValueError(f'Panel grid too small: rows*cols={nrows*ncols} < tracks={n_tracks}')

    if figsize is None:
        figsize = (ncols * cell_w, nrows * cell_h)

    def _draw_one(fig, axes, annotated: bool):
        if nrows == 1 and ncols == 1:
            axes_grid = np.array([[axes]])
        elif nrows == 1:
            axes_grid = np.array([axes])
        else:
            axes_grid = np.array(axes)

        for i, (tname, mat) in enumerate(track_mats_z):
            r = i // ncols
            c = i % ncols
            ax = axes_grid[r][c]

            sub = mat.loc[order]
            sns.heatmap(sub.values, cmap=cmap, yticklabels=False, xticklabels=False,
                        vmin=vmin, vmax=vmax, cbar=False, ax=ax)

            if annotated:
                ax.set_title(tname, fontsize=10)
                ax.set_xticks(region_spec.xticks_pos)
                ax.set_xticklabels(region_spec.xticks_labels, fontsize=8, rotation=0)

                for xp in [0, region_spec.bins_total]:
                    ax.plot([xp, xp], [0, len(order)], color='black', linewidth=0.6)
                for yp in [0, len(order)]:
                    ax.plot([0, region_spec.bins_total], [yp, yp], color='black', linewidth=0.6)
                for xp in region_spec.vlines_pos:
                    ax.plot([xp, xp], [0, len(order)], color='black', linewidth=0.6)
            else:
                ax.set_title('')
                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)

        for j in range(n_tracks, nrows * ncols):
            r = j // ncols
            c = j % ncols
            axes_grid[r][c].axis('off')

    # annotated
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    _draw_one(fig, axes, annotated=True)
    plt.tight_layout()
    plt.savefig(out_png_annot, dpi=dpi)
    plt.close(fig)

    # clean
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    _draw_one(fig, axes, annotated=False)
    plt.tight_layout()
    plt.savefig(out_png_clean, dpi=dpi)
    plt.close(fig)


def plot_single_heatmaps_png(track_mats_z: List[Tuple[str, pd.DataFrame]],
                             labels: pd.Series,
                             region_spec: RegionSpec,
                             out_dir: str,
                             cmap: str,
                             vmin: float,
                             vmax: float,
                             fig_x: float = 5.0,
                             fig_y: float = 7.0,
                             dpi: int = 200):
    """Save each track as an individual heatmap PNG (annot + clean).

    Output structure:
      out_dir/
        annot/<track>.png
        clean/<track>.png
    """

    order = labels.sort_values().index
    annot_dir = Path(out_dir) / 'annot'
    clean_dir = Path(out_dir) / 'clean'
    annot_dir.mkdir(parents=True, exist_ok=True)
    clean_dir.mkdir(parents=True, exist_ok=True)

    for tname, mat in track_mats_z:
        safe = _sanitize_filename(tname)
        sub = mat.loc[order]

        # annotated
        fig, ax = plt.subplots(1, 1, figsize=(fig_x, fig_y))
        sns.heatmap(sub.values, cmap=cmap, yticklabels=False, xticklabels=False,
                    vmin=vmin, vmax=vmax, cbar=False, ax=ax)
        ax.set_title(tname, fontsize=10)
        ax.set_xticks(region_spec.xticks_pos)
        ax.set_xticklabels(region_spec.xticks_labels, fontsize=8, rotation=0)
        for xp in [0, region_spec.bins_total]:
            ax.plot([xp, xp], [0, len(order)], color='black', linewidth=0.6)
        for yp in [0, len(order)]:
            ax.plot([0, region_spec.bins_total], [yp, yp], color='black', linewidth=0.6)
        for xp in region_spec.vlines_pos:
            ax.plot([xp, xp], [0, len(order)], color='black', linewidth=0.6)
        plt.tight_layout()
        fig.savefig(annot_dir / f'{safe}.png', dpi=dpi)
        plt.close(fig)

        # clean
        fig, ax = plt.subplots(1, 1, figsize=(fig_x, fig_y))
        sns.heatmap(sub.values, cmap=cmap, yticklabels=False, xticklabels=False,
                    vmin=vmin, vmax=vmax, cbar=False, ax=ax)
        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        plt.tight_layout()
        fig.savefig(clean_dir / f'{safe}.png', dpi=dpi)
        plt.close(fig)

def plot_cluster_expr_boxplot(expr_means: pd.DataFrame,
                              labels: pd.Series,
                              out_pdf: str,
                              ylim: Optional[Tuple[float, float]] = (0, 120)):
    """Expression boxplot by cluster."""
    m = expr_means.join(labels.rename('label'), how='inner')
    long = m.set_index('label').stack().reset_index()
    long.columns = ['label', 'treatment', 'FPKM']

    plt.figure(figsize=(8, 6), layout='constrained')
    sns.boxplot(x='label', y='FPKM', hue='treatment', data=long, palette='Set2', fliersize=False)
    if ylim is not None:
        plt.ylim(*ylim)
    sns.despine(trim=True)
    plt.legend(loc=[1.02, 0.35], frameon=False)
    plt.savefig(out_pdf)
    plt.close()


def plot_cluster_profiles(raw_mats: Dict[str, Dict[str, pd.DataFrame]],
                          marks: List[str],
                          treatments: List[str],
                          labels: pd.Series,
                          region_spec: RegionSpec,
                          out_pdf: str,
                          k: int,
                          profile_cols: int = 6):
    """Line profiles: for each cluster (row) x mark (col chunk), hue=treatment."""

    label_sorted = labels.sort_values()

    with PdfPages(out_pdf) as pdf:
        for start in range(0, len(marks), profile_cols):
            chunk = marks[start:start + profile_cols]
            ncols = len(chunk)
            fig, axes = plt.subplots(k, ncols, figsize=(ncols * 4.2, k * 2.8), sharex=True)
            if k == 1:
                axes = np.array([axes])
            if ncols == 1:
                axes = axes.reshape(k, 1)

            for ci, mark in enumerate(chunk):
                for clu in range(k):
                    ax = axes[clu][ci]
                    genes = label_sorted[label_sorted == clu].index

                    for tr in treatments:
                        mat = raw_mats[mark][tr]
                        g = genes.intersection(mat.index)
                        if len(g) == 0:
                            continue
                        y = mat.loc[g].mean(axis=0).values
                        ax.plot(np.arange(len(y)), y, label=tr)

                    for xp in region_spec.vlines_pos:
                        ax.axvline(xp, linestyle='--', color='grey', linewidth=1)

                    ax.set_xticks(region_spec.xticks_pos)
                    ax.set_xticklabels(region_spec.xticks_labels, fontsize=8, rotation=0)

                    if clu == 0:
                        ax.set_title(mark)
                    if ci == 0:
                        ax.set_ylabel(f"cluster {clu}")
                    else:
                        ax.set_ylabel('')

                    if (clu == 0) and (ci == ncols - 1):
                        ax.legend(loc=[1.02, 0.35], frameon=False)
                    else:
                        lg = ax.get_legend()
                        if lg:
                            lg.remove()

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


# -----------------------------
# Main pipeline
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            'Cluster genes using zscore-normalized ChIP/ATAC matrices (OmicsCanvas Step2 outputs) and plot:\n'
            '  (1) cluster heatmaps (panel + per-track),\n'
            '  (2) per-cluster mean profiles,\n'
            '  (3) per-cluster expression boxplots.'
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    ap.add_argument(
        '--matrix-dir',
        required=True,
        help=(
            'Directory containing Step2 matrices.\n'
            'Files are located by <matrix-dir>/<prefix><suffix>.\n'
            'Example: caculate_matrix'
        ),
    )

    ap.add_argument(
        '--in-group',
        required=True,
        help=(
            'Matrix file prefixes grouped by mark and treatment.\n'
            ' - Separate marks with ;\n'
            ' - Within a mark, separate treatments with ,\n'
            'Example:\n'
            '  G_H3K4me3,B_H3K4me3,R_H3K4me3;G_H3K27me3,B_H3K27me3,R_H3K27me3\n'
            'Each prefix resolves to: <matrix-dir>/<prefix><suffix>.'
        ),
    )
    ap.add_argument(
        '--in-ylabels',
        required=True,
        help=(
            'Mark names (one per ; block in --in-group).\n'
            'Comma-separated, used as plot titles/labels.\n'
            'Example: H3K4me3,H3K27me3'
        ),
    )
    ap.add_argument(
        '--in-names',
        required=True,
        help=(
            'Treatment names aligned to each prefix list in --in-group.\n'
            'Separate marks with ; and treatments with , (same structure as --in-group).\n'
            'Example: G,B,R;G,B,R'
        ),
    )

    ap.add_argument('--cluster-region', default='gene', choices=['TSS', 'gene', 'TES'],
                    help='Region used to build features and run KMeans.')
    ap.add_argument('--plot-regions', default='TSS,gene,TES',
                    help='Comma-separated regions to output plots for (e.g. "gene" or "TSS,gene,TES").')

    ap.add_argument('--naming', default='standard', choices=['standard', 'legacy'],
                    help='Matrix filename convention. standard uses *_tss_matrix.tsv / *_gene_profile_matrix.tsv / *_tes_matrix.tsv; '
                         'legacy uses *_start_matrix.txt / *_gene_body_matrix.txt / *_end_matrix.txt.')

    ap.add_argument('--suffix-tss', default=None,
                    help='Override TSS suffix (advanced). If omitted, inferred from --naming.')
    ap.add_argument('--suffix-gene', default=None,
                    help='Override gene suffix (advanced). If omitted, inferred from --naming.')
    ap.add_argument('--suffix-tes', default=None,
                    help='Override TES suffix (advanced). If omitted, inferred from --naming.')

    ap.add_argument('--distance', type=int, default=2000)
    ap.add_argument('--bins-start', type=int, default=100)
    ap.add_argument('--bins-end', type=int, default=100)
    ap.add_argument('--bins-gene-st', type=int, default=50)
    ap.add_argument('--bins-gene-body', type=int, default=100)
    ap.add_argument('--bins-gene-en', type=int, default=50)

    ap.add_argument('--scale-factor', type=float, default=10.0)
    ap.add_argument('--isoform-suffix', default='.1', help='Set empty to disable filtering.')

    ap.add_argument('--include-marks', default='', help='Comma-separated marks to include. Empty means all.')
    ap.add_argument('--exclude-marks', default='Input,Pol_II_Ser2P,Pol_II_Ser5P',
                    help='Comma-separated marks to exclude (from clustering + plots).')

    ap.add_argument('--k', type=int, default=8)
    ap.add_argument('--random-state', type=int, default=0)

    ap.add_argument('--cmap', default='RdBu_r')
    ap.add_argument('--vmin', type=float, default=-2)
    ap.add_argument('--vmax', type=float, default=2)

    ap.add_argument('--cols-per-row', type=int, default=9,
                    help='(compat) Panel columns per row. Prefer --panel-cols/--panel-rows for explicit grid.')
    ap.add_argument('--panel-cols', type=int, default=None,
                    help='Panel heatmap number of columns. Overrides --cols-per-row if set.')
    ap.add_argument('--panel-rows', type=int, default=None,
                    help='Panel heatmap number of rows. If set and --panel-cols not set, columns will be inferred.')
    ap.add_argument('--panel-cell-width', type=float, default=3.5,
                    help='Panel cell width (inches) when --panel-figsize is not provided.')
    ap.add_argument('--panel-cell-height', type=float, default=4.0,
                    help='Panel cell height (inches) when --panel-figsize is not provided.')
    ap.add_argument('--panel-figsize', default='',
                    help='Explicit panel figsize in inches, e.g. "30,15". Overrides cell width/height.')

    ap.add_argument('--profile-cols', type=int, default=6,
                    help='How many marks per page in the cluster profile PDF.')

    ap.add_argument('--heatmap-dpi', type=int, default=200,
                    help='DPI for PNG heatmaps (panel + per-track).')

    ap.add_argument('--save-single', dest='save_single', action='store_true', default=True,
                    help='Save each track heatmap as separate PNG files (annot+clean). Default: enabled.')
    ap.add_argument('--no-single', dest='save_single', action='store_false',
                    help='Disable per-track single heatmap outputs.')
    ap.add_argument('--single-fig-x', type=float, default=5.0,
                    help='Single-track heatmap width (inches).')
    ap.add_argument('--single-fig-y', type=float, default=7.0,
                    help='Single-track heatmap height (inches).')

    # ID mapping
    ap.add_argument('--id-mode', default='auto',
                    choices=['auto', 'none', 'drop_isoform', 'tx_to_gene', 'regex'],
                    help='How to normalize matrix IDs to match FPKM IDs. auto tries none/drop_isoform/tx_to_gene.')
    ap.add_argument('--id-regex', default='', help='Only used when --id-mode regex. Python regex pattern.')
    ap.add_argument('--id-repl', default='', help='Only used when --id-mode regex. Replacement string.')
    ap.add_argument('--collapse-duplicates', default='mean', choices=['mean', 'sum', 'first'],
                    help='If ID normalization creates duplicates, how to collapse rows.')

    ap.add_argument(
        '--fpkm',
        required=True,
        help=(
            "Expression table (TSV). The first column is the gene ID index; "
            "all remaining columns are numeric expression values (replicates)."
        ),
    )
    ap.add_argument('--expr-cols', required=True,
                    help='e.g. "G=0,1,2;B=3,4,5;R=6,7,8" (0-based indices)')
    ap.add_argument('--expr-ylim', default='0,120', help='e.g. "0,120". Empty to disable.')

    ap.add_argument('--out-prefix', default='histone_cluster')
    ap.add_argument('--outdir', default='histone_cluster_out')

    args = ap.parse_args()

    matrix_dir = args.matrix_dir
    if not matrix_dir.endswith('/'):
        matrix_dir += '/'

    os.makedirs(args.outdir, exist_ok=True)
    fig_dir = os.path.join(args.outdir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    # parse config strings
    marks_all, prefixes_per_mark, treatments_per_mark = parse_in_group(args.in_group, args.in_ylabels, args.in_names)

    include_marks = set(_split_nonempty(args.include_marks, ',')) if args.include_marks.strip() else set(marks_all)
    exclude_marks = set(_split_nonempty(args.exclude_marks, ','))

    marks = [m for m in marks_all if (m in include_marks) and (m not in exclude_marks)]
    if not marks:
        raise ValueError('No marks selected after include/exclude filtering.')

    # assume consistent treatments across marks (same order)
    treatments = treatments_per_mark[0]

    defaults = {
        'standard': {'TSS': '_tss_matrix.tsv', 'gene': '_gene_profile_matrix.tsv', 'TES': '_tes_matrix.tsv'},
        'legacy':   {'TSS': '_start_matrix.txt', 'gene': '_gene_body_matrix.txt', 'TES': '_end_matrix.txt'},
    }
    suffix_map = defaults[args.naming].copy()
    if args.suffix_tss is not None: suffix_map['TSS'] = args.suffix_tss
    if args.suffix_gene is not None: suffix_map['gene'] = args.suffix_gene
    if args.suffix_tes is not None: suffix_map['TES'] = args.suffix_tes

    plot_regions = _split_nonempty(args.plot_regions, ',')
    for r in plot_regions:
        if r not in suffix_map:
            raise ValueError(f"Bad region in --plot-regions: {r}")

    regions_needed = sorted(set([args.cluster_region] + plot_regions))

    # region specs
    region_specs: Dict[str, RegionSpec] = {}
    for r in regions_needed:
        region_specs[r] = build_region_spec(
            region=r,
            distance=args.distance,
            bins_start=args.bins_start,
            bins_end=args.bins_end,
            bins_gene_st=args.bins_gene_st,
            bins_gene_body=args.bins_gene_body,
            bins_gene_en=args.bins_gene_en,
        )

    # expression
    expr = pd.read_csv(args.fpkm, sep='\t', index_col=0)
    expr.index = expr.index.astype(str)
    expr_cols_map = parse_expr_cols(args.expr_cols)

    expr_means_full = pd.DataFrame(index=expr.index)
    for tr, cols in expr_cols_map.items():
        if max(cols) >= expr.shape[1]:
            raise ValueError(f"expr cols out of range for treatment {tr}: {cols}")
        expr_means_full[tr] = pd.to_numeric(expr.iloc[:, cols].mean(axis=1), errors='coerce').fillna(0.0)

    # ID mapping (decide once, then apply to ALL matrices)
    id_mode = args.id_mode

    # load matrices: raw_by_region[region][mark][treatment]
    raw_by_region: Dict[str, Dict[str, Dict[str, pd.DataFrame]]] = {r: {} for r in regions_needed}

    mark_to_i = {m: i for i, m in enumerate(marks_all)}

    # first pass: read one matrix to auto-pick id-mode (if auto)
    probe_matrix_ids: List[str] = []
    if id_mode == 'auto':
        # pick the first selected mark + first treatment + cluster_region
        probe_mark = marks[0]
        i = mark_to_i[probe_mark]
        probe_prefix = prefixes_per_mark[i][0]
        probe_suffix = suffix_map[args.cluster_region]
        probe_path = os.path.join(matrix_dir, probe_prefix + probe_suffix)
        if not os.path.exists(probe_path):
            raise FileNotFoundError(f"Missing matrix (probe for id-mode): {probe_path}")
        probe_df = read_matrix(probe_path, scale_factor=args.scale_factor, isoform_suffix=args.isoform_suffix)
        probe_matrix_ids = probe_df.index.astype(str).tolist()
        id_mode = pick_best_id_mode(probe_matrix_ids, expr_means_full.index.tolist(), allow=['none', 'drop_isoform', 'tx_to_gene'])
        print(f"[INFO] --id-mode auto selected: {id_mode}")

    id_func = _make_id_transform(id_mode, args.id_regex, args.id_repl)

    # load all matrices and apply ID transform + collapse
    for r in regions_needed:
        suffix = suffix_map[r]
        spec = region_specs[r]

        for mark in marks:
            i = mark_to_i[mark]
            raw_by_region[r][mark] = {}
            for j, tr in enumerate(treatments_per_mark[i]):
                prefix = prefixes_per_mark[i][j]
                path = os.path.join(matrix_dir, prefix + suffix)
                if not os.path.exists(path):
                    raise FileNotFoundError(f"Missing matrix: {path}")

                df = read_matrix(path, scale_factor=args.scale_factor, isoform_suffix=args.isoform_suffix)
                df = apply_id_transform_and_collapse(df, id_func=id_func, collapse=args.collapse_duplicates)

                if df.shape[1] != spec.bins_total:
                    print(f"[WARN] {r} {mark}|{tr} bins={df.shape[1]} != expected {spec.bins_total}. Continue.")
                raw_by_region[r][mark][tr] = df

    # strict gene intersection across all required regions (after ID normalization)
    all_dfs = []
    for r in regions_needed:
        for mark in marks:
            for tr in treatments:
                all_dfs.append(raw_by_region[r][mark][tr])

    common_genes_all = gene_intersection(all_dfs)
    if not common_genes_all:
        raise ValueError('No common genes across selected tracks/regions (after ID normalization + isoform filtering).')

    # build features from cluster-region only
    cluster_raw = raw_by_region[args.cluster_region]

    features = pd.DataFrame(index=common_genes_all)

    for mark in marks:
        for tr in treatments:
            df = cluster_raw[mark][tr].loc[common_genes_all]
            zdf = zscore_by_gene(df)
            zdf2 = zdf.copy()
            zdf2.columns = [f"{mark}|{tr}|{c}" for c in zdf2.columns]
            features = pd.concat([features, zdf2], axis=1)

    features = features.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how='any')
    genes_used = features.index.tolist()
    if len(genes_used) == 0:
        raise ValueError('All genes dropped after NaN filtering; check matrices.')

    # clustering
    km = KMeans(n_clusters=args.k, random_state=args.random_state, n_init=10)
    labels = pd.Series(km.fit_predict(features.values), index=features.index, name='label')

    labels_out = os.path.join(args.outdir, f"{args.out_prefix}_cluster_region-{args.cluster_region}_k{args.k}_labels.tsv")
    labels.sort_values().to_csv(labels_out, sep='\t', header=True)

    feat_out = os.path.join(args.outdir, f"{args.out_prefix}_feature_matrix_{args.cluster_region}_zscore.tsv.gz")
    with gzip.open(feat_out, 'wt', encoding='utf-8') as w:
        features.to_csv(w, sep='\t')

    # expression means (subset to clustered genes)
    expr_means = expr_means_full.loc[expr_means_full.index.intersection(genes_used)]

    if expr_means.shape[0] == 0:
        # show a helpful debug message
        ex_m = probe_matrix_ids[:5] if probe_matrix_ids else []
        ex_m = [id_func(x) for x in ex_m]
        ex_e = expr_means_full.index[:5].tolist()
        raise ValueError(
            "Expression join is empty (no overlapping IDs between matrices and FPKM).\n"
            f"  id_mode used: {id_mode}\n"
            f"  example matrix IDs (mapped): {ex_m}\n"
            f"  example expr IDs: {ex_e}\n"
            "Fix by using --id-mode (none/drop_isoform/tx_to_gene/regex) to match your ID systems."
        )

    ylim = None
    if args.expr_ylim.strip():
        a, b = args.expr_ylim.split(',')
        ylim = (float(a), float(b))

    expr_pdf = os.path.join(fig_dir, f"{args.out_prefix}_cluster_expr_boxplot_k{args.k}.pdf")
    plot_cluster_expr_boxplot(expr_means, labels, expr_pdf, ylim=ylim)

    # plots for each requested region
    for r in plot_regions:
        spec = region_specs[r]
        raw_r = raw_by_region[r]

        # zscore mats for heatmaps in region r (subset to genes_used)
        track_mats_z: List[Tuple[str, pd.DataFrame]] = []
        for mark in marks:
            for tr in treatments:
                df = raw_r[mark][tr].loc[genes_used]
                zdf = zscore_by_gene(df)
                track_mats_z.append((f"{mark}|{tr}", zdf))

        heat_png_annot = os.path.join(fig_dir, f"{args.out_prefix}_{r}_panel_heatmaps_annot_k{args.k}.png")
        heat_png_clean = os.path.join(fig_dir, f"{args.out_prefix}_{r}_panel_heatmaps_clean_k{args.k}.png")

        panel_figsize = _parse_figsize_str(args.panel_figsize)
        panel_cols = args.panel_cols if args.panel_cols is not None else args.cols_per_row

        plot_panel_heatmaps_png(
            track_mats_z=track_mats_z,
            labels=labels,
            region_spec=spec,
            out_png_annot=heat_png_annot,
            out_png_clean=heat_png_clean,
            cmap=args.cmap,
            vmin=args.vmin,
            vmax=args.vmax,
            panel_cols=panel_cols,
            panel_rows=args.panel_rows,
            cell_w=args.panel_cell_width,
            cell_h=args.panel_cell_height,
            figsize=panel_figsize,
            dpi=args.heatmap_dpi,
        )

        # per-track outputs (into a folder)
        if args.save_single:
            single_dir = os.path.join(fig_dir, 'single_tracks', r)
            plot_single_heatmaps_png(
                track_mats_z=track_mats_z,
                labels=labels,
                region_spec=spec,
                out_dir=single_dir,
                cmap=args.cmap,
                vmin=args.vmin,
                vmax=args.vmax,
                fig_x=args.single_fig_x,
                fig_y=args.single_fig_y,
                dpi=args.heatmap_dpi,
            )

        # profiles (compact -> PDF ok)
        raw_used: Dict[str, Dict[str, pd.DataFrame]] = {}
        for mark in marks:
            raw_used[mark] = {}
            for tr in treatments:
                raw_used[mark][tr] = raw_r[mark][tr].loc[genes_used]

        prof_pdf = os.path.join(fig_dir, f"{args.out_prefix}_{r}_cluster_profiles_k{args.k}.pdf")
        plot_cluster_profiles(
            raw_mats=raw_used,
            marks=marks,
            treatments=treatments,
            labels=labels,
            region_spec=spec,
            out_pdf=prof_pdf,
            k=args.k,
            profile_cols=args.profile_cols,
        )

    print("[OK] Done")
    print(f"  labels:   {labels_out}")
    print(f"  feature:  {feat_out}")
    print(f"  expr:     {expr_pdf}")
    for r in plot_regions:
        print(f"  heatmap({r}) annotated: {os.path.join(fig_dir, f'{args.out_prefix}_{r}_panel_heatmaps_annot_k{args.k}.png')}")
        print(f"  heatmap({r}) clean:     {os.path.join(fig_dir, f'{args.out_prefix}_{r}_panel_heatmaps_clean_k{args.k}.png')}")
        print(f"  profile({r}):           {os.path.join(fig_dir, f'{args.out_prefix}_{r}_cluster_profiles_k{args.k}.pdf')}")


if __name__ == '__main__':
    main()
