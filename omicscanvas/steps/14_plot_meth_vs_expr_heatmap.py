#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
meth_vs_expr_heatmap.py
=======================

根据表达量（FPKM/TPM/Counts）对基因排序分箱，然后计算每个表达分箱在基因结构（start/body/end）
上的甲基化比例（sum(me)/sum(al)）并绘制热图。

输入甲基化文件（与你当前一致）：
  <meth_dir>/<sample>_<CX>_gene_start_300.txt
  <meth_dir>/<sample>_<CX>_gene_gene_300.txt
  <meth_dir>/<sample>_<CX>_gene_end_300.txt

甲基化表要求至少含列：
  name   po   me   al
说明：
  - 你读文件时常用 index_col=0，会导致 name 在 index 里；本脚本自动兼容。
  - po 可能不是从 0/1 开始，比如 1..100 或其它起点；本脚本会根据每段 po_min/po_max
    自动归一化并把 gene/end 段累加偏移到全局坐标。

表达矩阵：
  - tsv，index 为 gene ID
  - 使用 --fpkm-cols 指定参与均值的列（从 0 开始，不含 index）

输出：
  <out-prefix>_<sample>_<CX>_meth_vs_expr_heatmap.pdf

示例（按你刚才的命令）：
  python meth_vs_expr_heatmap.py \
    --sample LSH10 \
    --cx CG \
    --meth-dir gene \
    --fpkm FPKM.txt \
    --fpkm-cols 3,4,5 \
    --none-bins 10 \
    --exp-bins 90 \
    --distance 2000 \
    --suffix-start _gene_start_300.txt \
    --suffix-gene  _gene_gene_300.txt \
    --suffix-end   _gene_end_300.txt \
    --scale-mode ratio \
    --range 1 \
    --cmap RdBu_r \
    --out-prefix LSH10
"""

import os
import argparse
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# -----------------------------
# Basic utils
# -----------------------------

def ensure_dir(d: str) -> str:
    if d and (not d.endswith("/")):
        d += "/"
    return d


def parse_int_list(s: str) -> List[int]:
    parts = [x.strip() for x in s.split(",") if x.strip()]
    return [int(x) for x in parts]


def parse_float_pair(s: str, name: str) -> Tuple[float, float]:
    a, b = [x.strip() for x in s.split(",")]
    return float(a), float(b)


# -----------------------------
# Read methylation tables
# -----------------------------

def read_meth_table(path: str) -> pd.DataFrame:
    """
    读取甲基化表（tsv）
    兼容两种常见格式：
      A) 文件有列：name po me al（name 在列里）
      B) 你用 index_col=0 时：name 在 index 里，列为 po me al

    返回列固定为：name, po, me, al
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")

    # 先按你常用的方式读（index_col=0），兼容 B
    df = pd.read_csv(path, sep="\t", index_col=0)

    # 若 name 不在列里，则把 index 当作 name（兼容 B）
    if "name" not in df.columns:
        df = df.copy()
        df["name"] = df.index.astype(str)

    required = {"name", "po", "me", "al"}
    miss = required - set(df.columns)
    if miss:
        raise ValueError(f"{path} missing columns: {sorted(miss)}; got columns={list(df.columns)}")

    out = df[["name", "po", "me", "al"]].copy()
    out["name"] = out["name"].astype(str)
    out["po"] = pd.to_numeric(out["po"], errors="coerce")
    out["me"] = pd.to_numeric(out["me"], errors="coerce")
    out["al"] = pd.to_numeric(out["al"], errors="coerce")
    out = out.dropna(subset=["name", "po", "me", "al"])
    out["po"] = out["po"].astype(int)

    return out


def infer_po_range(df: pd.DataFrame) -> Tuple[int, int, int]:
    """
    返回 (po_min, po_max, seg_bins) 其中 seg_bins = po_max - po_min + 1
    """
    po_vals = pd.to_numeric(df["po"], errors="coerce").dropna()
    if po_vals.empty:
        raise ValueError("po is not numeric / empty.")
    po_min = int(po_vals.min())
    po_max = int(po_vals.max())
    seg_bins = po_max - po_min + 1
    if seg_bins <= 0:
        raise ValueError(f"Invalid seg_bins={seg_bins} from po_min={po_min}, po_max={po_max}")
    return po_min, po_max, seg_bins


def normalize_and_offset_po(df: pd.DataFrame, po_min: int, offset: int) -> pd.DataFrame:
    """
    将段内 po 归一化到 0..(seg_bins-1)，然后加 offset 变成全局坐标
      po_global = (po - po_min) + offset
    """
    d = df.copy()
    d["po"] = (d["po"].astype(int) - int(po_min)) + int(offset)
    return d


def load_and_merge_methylation(meth_dir: str,
                               sample: str,
                               cx: str,
                               suffix_start: str,
                               suffix_gene: str,
                               suffix_end: str,
                               allow_unequal_segments: bool = False) -> Tuple[pd.DataFrame, List[int], List[int]]:
    """
    读取三段并拼接成全局 po。
    返回：
      merged_df: 包含 name, po(全局), me, al
      seg_bins_list: [seg1, seg2, seg3]
      boundaries: [b1, b2, total] 其中 b1=seg1, b2=seg1+seg2, total=seg1+seg2+seg3
        这些用于画竖线和刻度
    """
    meth_dir = ensure_dir(meth_dir)
    f1 = f"{meth_dir}{sample}_{cx}{suffix_start}"
    f2 = f"{meth_dir}{sample}_{cx}{suffix_gene}"
    f3 = f"{meth_dir}{sample}_{cx}{suffix_end}"

    on1 = read_meth_table(f1)
    on2 = read_meth_table(f2)
    on3 = read_meth_table(f3)

    po1_min, po1_max, seg1 = infer_po_range(on1)
    po2_min, po2_max, seg2 = infer_po_range(on2)
    po3_min, po3_max, seg3 = infer_po_range(on3)

    if not (seg1 == seg2 == seg3) and (not allow_unequal_segments):
        raise ValueError(
            "Segment bins are not equal.\n"
            f" start: seg={seg1}, po=[{po1_min},{po1_max}]\n"
            f" gene : seg={seg2}, po=[{po2_min},{po2_max}]\n"
            f" end  : seg={seg3}, po=[{po3_min},{po3_max}]\n"
            "If you still want to proceed, add --allow-unequal-segments"
        )

    # 累加偏移（允许不等长时也能拼接，只是 TSS/TES 边界会跟着变）
    b1 = seg1
    b2 = seg1 + seg2
    total = seg1 + seg2 + seg3

    on1n = normalize_and_offset_po(on1, po1_min, 0)
    on2n = normalize_and_offset_po(on2, po2_min, b1)
    on3n = normalize_and_offset_po(on3, po3_min, b2)

    merged = pd.concat([on1n, on2n, on3n], ignore_index=True)

    # 预聚合：减少重复行（你的数据每个 gene/po 会多行）
    merged = merged.groupby(["name", "po"], as_index=False)[["me", "al"]].sum()

    print(f"[INFO] start po=[{po1_min},{po1_max}] seg={seg1} -> global [0,{seg1-1}]")
    print(f"[INFO] gene  po=[{po2_min},{po2_max}] seg={seg2} -> global [{b1},{b2-1}]")
    print(f"[INFO] end   po=[{po3_min},{po3_max}] seg={seg3} -> global [{b2},{total-1}]")
    print(f"[INFO] total_bins={total}")

    return merged, [seg1, seg2, seg3], [b1, b2, total]


# -----------------------------
# Expression loading & binning
# -----------------------------

def load_expression(expr_path: str, cols: List[int]) -> pd.Series:
    """
    读取表达矩阵并对指定列求均值（index 为 gene）
    """
    if not os.path.exists(expr_path):
        raise FileNotFoundError(f"Not found: {expr_path}")

    df = pd.read_csv(expr_path, sep="\t", index_col=0)
    if df.shape[1] == 0:
        raise ValueError("Expression table has no columns.")

    max_col = df.shape[1] - 1
    for c in cols:
        if c < 0 or c > max_col:
            raise ValueError(f"fpkm-col {c} out of range (0..{max_col})")

    expr = df.iloc[:, cols].mean(axis=1)
    expr = pd.to_numeric(expr, errors="coerce").fillna(0.0)
    expr.index = expr.index.astype(str)
    return expr


def bin_genes(expr: pd.Series, none_bins: int, exp_bins: int, zero_threshold: float = 0.0) -> List[List[str]]:
    """
    将 overlap 基因按表达量排序分箱：
      - expr <= threshold -> none_bins 份
      - expr  > threshold -> exp_bins 份
    用 np.array_split：不会丢尾巴
    """
    expr_sorted = expr.sort_values(ascending=True)
    none_part = expr_sorted[expr_sorted <= zero_threshold]
    exp_part  = expr_sorted[expr_sorted >  zero_threshold]

    groups = []

    if none_bins > 0:
        none_splits = np.array_split(none_part.index.to_numpy(), none_bins) if len(none_part) > 0 else [np.array([]) for _ in range(none_bins)]
        for arr in none_splits:
            groups.append(arr.tolist())

    if exp_bins > 0:
        exp_splits = np.array_split(exp_part.index.to_numpy(), exp_bins) if len(exp_part) > 0 else [np.array([]) for _ in range(exp_bins)]
        for arr in exp_splits:
            groups.append(arr.tolist())

    return groups


def build_heat_matrix(meth_df: pd.DataFrame, groups: List[List[str]], total_bins: int, fillna_value: Optional[float]) -> np.ndarray:
    """
    对每个表达分箱 group 计算每个全局 po 的 ratio = sum(me)/sum(al)
    """
    full_index = pd.Index(range(total_bins), name="po")
    rows = []

    for genes in groups:
        if len(genes) == 0:
            row = np.full(total_bins, np.nan)
            rows.append(row)
            continue

        sub = meth_df[meth_df["name"].isin(genes)]
        if sub.shape[0] == 0:
            row = np.full(total_bins, np.nan)
            rows.append(row)
            continue

        agg = sub.groupby("po")[["me", "al"]].sum().reindex(full_index, fill_value=0.0)
        me = agg["me"].to_numpy()
        al = agg["al"].to_numpy()

        ratio = np.divide(me, al, out=np.full(total_bins, np.nan), where=(al != 0))

        if fillna_value is not None:
            ratio = np.where(np.isfinite(ratio), ratio, fillna_value)

        rows.append(ratio)

    return np.vstack(rows)


# -----------------------------
# vmin/vmax & cmap
# -----------------------------

def default_cmap(scale_mode: str) -> str:
    # ratio/quantile：更常见 sequential；diverging：RdBu_r
    if scale_mode == "diverging":
        return "RdBu_r"
    return "viridis"


def compute_vmin_vmax(heat: np.ndarray,
                      scale_mode: str,
                      user_range: Optional[float],
                      user_vmin: Optional[float],
                      user_vmax: Optional[float],
                      q_low: float,
                      q_high: float,
                      clamp_ratio_to_1: bool = True) -> Tuple[float, float]:
    """
    优先级：
      1) 同时给 vmin+vmax -> 直接用
      2) 给 range -> ratio/quantile: [0, range]; diverging: [-range, +range]
      3) quantile -> 用分位数
      4) ratio -> vmin=0, vmax=99%分位（可 clamp<=1）
      5) diverging -> ±abs的99%分位
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
            vmax =  r if user_vmax is None else float(user_vmax)
        else:
            vmin = 0.0 if user_vmin is None else float(user_vmin)
            vmax = r   if user_vmax is None else float(user_vmax)
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

    if scale_mode == "ratio":
        vmin = 0.0 if user_vmin is None else float(user_vmin)
        vmax = float(np.nanquantile(data, 0.99)) if user_vmax is None else float(user_vmax)
        if clamp_ratio_to_1:
            vmax = min(vmax, 1.0)
        if vmax <= vmin:
            vmax = vmin + 1e-6
        return vmin, vmax

    if scale_mode == "diverging":
        peak = float(np.nanquantile(np.abs(data), 0.99))
        vmin = -peak if user_vmin is None else float(user_vmin)
        vmax =  peak if user_vmax is None else float(user_vmax)
        if vmax <= vmin:
            vmax = vmin + 1e-6
        return vmin, vmax

    # fallback
    vmin = float(np.nanmin(data)) if user_vmin is None else float(user_vmin)
    vmax = float(np.nanmax(data)) if user_vmax is None else float(user_vmax)
    if vmax <= vmin:
        vmax = vmin + 1e-6
    return vmin, vmax


# -----------------------------
# Plot
# -----------------------------

def plot_heatmap(heat: np.ndarray,
                 out_pdf: str,
                 distance: int,
                 boundaries: List[int],
                 cmap: str,
                 vmin: float,
                 vmax: float,
                 show_box: bool,
                 ytick_step: Optional[int] = None):
    """
    boundaries: [b1, b2, total]
      - start: 0..b1-1
      - gene : b1..b2-1
      - end  : b2..total-1
    """
    nrows, total_bins = heat.shape
    b1, b2, total = boundaries

    # figure size: 行多就拉高
    fig_h = max(6, nrows / 10.0)
    fig_w = 8
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = plt.gca()

    sns.heatmap(
        heat,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        cbar=True
    )

    # x ticks：放在 cell center（+0.5）
    xticks = [0, b1 - 1, b2 - 1, total - 1]
    xlabels = [f"-{distance}bp", "TSS", "TES", f"+{distance}bp"]
    ax.set_xticks([x + 0.5 for x in xticks])
    ax.set_xticklabels(xlabels, rotation=0)

    # y ticks
    if ytick_step is None:
        ytick_step = max(1, nrows // 10)
    yticks = list(range(0, nrows, ytick_step))
    ax.set_yticks([y + 0.5 for y in yticks])
    ax.set_yticklabels(yticks, rotation=0)

    # 竖线：段边界（画在边界处）
    ax.axvline(b1, linestyle="--", color="black", linewidth=1)
    ax.axvline(b2, linestyle="--", color="black", linewidth=1)

    if show_box:
        ax.plot([0, total_bins], [0, 0], color="black", linewidth=1)
        ax.plot([0, total_bins], [nrows, nrows], color="black", linewidth=1)
        ax.plot([0, 0], [0, nrows], color="black", linewidth=1)
        ax.plot([total_bins, total_bins], [0, nrows], color="black", linewidth=1)

    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Methylation-vs-expression heatmap. Genes are ranked by expression, binned, then methylation ratio (sum(me)/sum(al)) is computed per bin along the merged profile. Input methylation files: <meth_dir>/<sample>_<cx><suffix-*>. Example: CX_gene/SRR9321764_CHH_upstream_bins50.tsv (use --suffix-start _upstream_bins50.tsv)."
    )

    ap.add_argument(
        "--sample", required=True,
        help=(
            "Sample ID used to build filenames. Example: SRR9321764. "
            "Files are searched as <meth_dir>/<sample>_<cx><suffix>."
        )
    )
    ap.add_argument(
        "--cx", required=True,
        help="Methylation context, e.g. CG / CHG / CHH."
    )
    ap.add_argument(
        "--meth-dir", default="gene",
        help="Directory containing methylation segment tables (default: gene/)."
    )

    ap.add_argument(
        "--suffix-start", default="_gene_start_300.txt",
        help=(
            "Suffix for upstream/TSS segment file. "
            "Example (your new naming): _upstream_bins50.tsv"
        )
    )
    ap.add_argument(
        "--suffix-gene", default="_gene_gene_300.txt",
        help=(
            "Suffix for gene body segment file. "
            "Example (your new naming): _body_bins100.tsv"
        )
    )
    ap.add_argument(
        "--suffix-end", default="_gene_end_300.txt",
        help=(
            "Suffix for downstream/TES segment file. "
            "Example (your new naming): _downstream_bins50.tsv"
        )
    )
    ap.add_argument(
        "--allow-unequal-segments", action="store_true",
        help=(
            "Allow different bin counts across start/gene/end segments. "
            "By default, the script requires seg_bins(start)=seg_bins(gene)=seg_bins(end)."
        )
    )

    ap.add_argument(
        "--fpkm", required=True,
        help=(
            "Expression table (TSV). Row index must be gene IDs. "
            "The script averages columns selected by --fpkm-cols."
        )
    )
    ap.add_argument(
        "--fpkm-cols", default="0,1,2",
        help=(
            "0-based column indices (excluding the index column) to average for expression. "
            "Example: 3,4,5"
        )
    )

    ap.add_argument("--none-bins", type=int, default=10, help="Number of bins for genes with expression <= threshold.")
    ap.add_argument("--exp-bins", type=int, default=90, help="Number of bins for genes with expression  > threshold.")
    ap.add_argument("--zero-threshold", type=float, default=0.0, help="Expression <= this value is treated as 'none' (default: 0).")

    ap.add_argument(
        "--distance", type=int, default=2000,
        help=(
            "Distance (bp) shown on x-axis labels (e.g. -2000bp / +2000bp). "
            "This does NOT affect po/bin merging."
        )
    )

    ap.add_argument(
        "--scale-mode", default="ratio", choices=["ratio", "diverging", "quantile"],
        help=(
            "How to set color scale. "
            "ratio: methylation ratio in [0,1]; diverging: symmetric around 0; quantile: based on quantiles."
        )
    )
    ap.add_argument(
        "--range", type=float, default=None,
        help=(
            "Preset range for color scale. For 'ratio', range=1 means [0,1]. "
            "For 'diverging', range=2 means [-2,2]. If you set --vmin/--vmax, this is ignored."
        )
    )
    ap.add_argument("--vmin", type=float, default=None, help="Manual vmin (optional).")
    ap.add_argument("--vmax", type=float, default=None, help="Manual vmax (optional).")
    ap.add_argument("--quantiles", default="0.01,0.99", help="Quantile bounds for scale-mode=quantile, e.g. 0.05,0.95")
    ap.add_argument("--no-clamp1", action="store_true", help="(ratio mode) Do not clamp vmax to <= 1.")

    ap.add_argument(
        "--cmap", default=None,
        help=(
            "Heatmap colormap name. If not provided, a default cmap is chosen based on --scale-mode."
        )
    )

    ap.add_argument(
        "--fillna", type=float, default=None,
        help=(
            "Fill value for missing bins (NaN). If not set, NaNs remain and may appear blank."
        )
    )
    ap.add_argument("--no-box", action="store_true", help="Do not draw the outer border.")
    ap.add_argument("--ytick-step", type=int, default=None, help="Y-axis tick step. If not set, the script auto-selects.")

    ap.add_argument(
        "--out-prefix", default=None,
        help="Output prefix. Default: use --sample."
    )

    args = ap.parse_args()

    # 1) 甲基化三段拼接（按 po 范围归一化 + 累加偏移）
    meth_df, segs, boundaries = load_and_merge_methylation(
        meth_dir=args.meth_dir,
        sample=args.sample,
        cx=args.cx,
        suffix_start=args.suffix_start,
        suffix_gene=args.suffix_gene,
        suffix_end=args.suffix_end,
        allow_unequal_segments=args.allow_unequal_segments
    )
    total_bins = boundaries[-1]

    # 2) 表达读取
    cols = parse_int_list(args.fpkm_cols)
    expr = load_expression(args.fpkm, cols)

    # 3) overlap genes
    meth_genes = set(meth_df["name"].unique())
    expr_genes = set(expr.index)
    overlap = sorted(meth_genes & expr_genes)
    if len(overlap) == 0:
        raise ValueError("No overlap genes between methylation and expression table.")
    expr_overlap = expr.loc[overlap]

    # 4) 表达分箱
    groups = bin_genes(expr_overlap, args.none_bins, args.exp_bins, args.zero_threshold)

    # 5) heat matrix
    heat = build_heat_matrix(meth_df, groups, total_bins, fillna_value=args.fillna)

    # 6) vmin/vmax/cmap
    q_low, q_high = parse_float_pair(args.quantiles, "--quantiles")
    vmin, vmax = compute_vmin_vmax(
        heat=heat,
        scale_mode=args.scale_mode,
        user_range=args.range,
        user_vmin=args.vmin,
        user_vmax=args.vmax,
        q_low=q_low,
        q_high=q_high,
        clamp_ratio_to_1=(not args.no_clamp1)
    )
    cmap = args.cmap if args.cmap is not None else default_cmap(args.scale_mode)

    # 7) 输出
    out_prefix = args.out_prefix if args.out_prefix else args.sample
    out_pdf = f"{out_prefix}_{args.sample}_{args.cx}_meth_vs_expr_heatmap.pdf"

    plot_heatmap(
        heat=heat,
        out_pdf=out_pdf,
        distance=args.distance,
        boundaries=boundaries,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        show_box=(not args.no_box),
        ytick_step=args.ytick_step
    )

    print(f"[OK] Saved: {out_pdf}")
    print(f"[INFO] overlap genes: {len(overlap)}")
    print(f"[INFO] heat shape: {heat.shape} (rows={args.none_bins+args.exp_bins}, cols={total_bins})")
    print(f"[INFO] scale-mode={args.scale_mode} cmap={cmap} vmin={vmin:.6g} vmax={vmax:.6g}")


if __name__ == "__main__":
    main()
