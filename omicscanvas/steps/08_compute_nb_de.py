#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
import pandas as pd

from statsmodels.discrete.discrete_model import NegativeBinomial
from statsmodels.stats.multitest import multipletests


def eprint(msg: str):
    print(msg, file=sys.stderr)


def read_counts(path: str) -> pd.DataFrame:
    """
    Expect genes x samples:
      first column: gene_id
      other columns: sample names
    Values must be raw counts (non-negative).
    """
    df = pd.read_csv(path, sep=None, engine="python")
    if df.shape[1] < 3:
        raise ValueError("Counts table should have >= 3 columns: gene_id + >=2 samples.")
    df = df.set_index(df.columns[0])
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    df = df.apply(pd.to_numeric, errors="coerce")
    if df.isna().any().any():
        bad = int(df.isna().sum().sum())
        raise ValueError(f"Counts contain NA/non-numeric entries: {bad}.")
    if (df < 0).any().any():
        raise ValueError("Counts contain negative values. Raw counts must be non-negative.")

    # DE expects integers; allow floats that are integers
    if not np.all(np.mod(df.values, 1) == 0):
        eprint("[WARN] Counts contain non-integers; rounding to nearest int. "
               "Make sure you are using raw counts, not TPM/FPKM.")
    df = df.round().astype(int)
    return df


def read_meta(path: str, sample_col: str, condition_col: str) -> pd.DataFrame:
    meta = pd.read_csv(path, sep=None, engine="python")
    if sample_col not in meta.columns:
        raise ValueError(f"Metadata missing sample column '{sample_col}'.")
    if condition_col not in meta.columns:
        raise ValueError(f"Metadata missing condition column '{condition_col}'.")
    meta[sample_col] = meta[sample_col].astype(str)
    meta[condition_col] = meta[condition_col].astype(str)
    meta = meta.set_index(sample_col)
    meta.index.name = "sample"
    return meta


def deseq_median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    """
    Median-of-ratios size factors (DESeq-style).
    - Compute per-gene geometric mean (robust: uses only positive counts).
    - For each sample, size factor = median(count_ij / geoMean_i) over genes with geoMean>0 and count_ij>0.
    """
    X = counts.to_numpy(dtype=float)  # genes x samples

    # geometric mean per gene using positive counts only
    with np.errstate(divide="ignore", invalid="ignore"):
        logX = np.where(X > 0, np.log(X), np.nan)
    gm = np.exp(np.nanmean(logX, axis=1))  # genes
    valid_gene = np.isfinite(gm) & (gm > 0)

    if valid_gene.sum() < 10:
        raise RuntimeError("Too few genes with positive counts to compute size factors.")

    # ratios per sample
    ratios = X[valid_gene, :] / gm[valid_gene, None]
    # take median over positive ratios only
    sf = []
    for j in range(ratios.shape[1]):
        r = ratios[:, j]
        r = r[np.isfinite(r) & (r > 0)]
        if len(r) == 0:
            sf.append(np.nan)
        else:
            sf.append(np.median(r))
    sf = np.array(sf, dtype=float)

    if not np.all(np.isfinite(sf)) or np.any(sf <= 0):
        raise RuntimeError("Failed to compute valid size factors. Check counts matrix.")

    # Normalize so geometric mean of size factors = 1 (optional but nice)
    sf = sf / np.exp(np.mean(np.log(sf)))
    return pd.Series(sf, index=counts.columns, name="size_factor")


def fit_nb_for_gene(y: np.ndarray, x: np.ndarray, offset: np.ndarray):
    """
    Fit NB regression for one gene:
      y ~ NB(mean=exp(offset + x*b), dispersion estimated)
    Using statsmodels.discrete NegativeBinomial (NB2) MLE.
    Returns: (beta_condition, se, pvalue, alpha)
    """
    try:
        model = NegativeBinomial(y, x, offset=offset)
        # disp=False prevents printing
        res = model.fit(disp=False, maxiter=100)
        # x columns: [intercept, condition]
        beta = float(res.params[1])
        se = float(res.bse[1])
        pval = float(res.pvalues[1])
        # statsmodels stores alpha as last parameter in some parameterizations;
        # but res.model is NB2; use res.model._dispersion or res.params? We'll try robust extraction:
        alpha = getattr(res, "params", None)
        # Safer: use res.model if available
        # statsmodels Discrete NB stores ln(alpha) or alpha depending on version;
        # We'll expose res.scale if alpha not easily accessible.
        # Here, use res.model if has attribute 'alpha' or 'lnalpha'.
        a = None
        if hasattr(res.model, "alpha"):
            a = float(res.model.alpha)
        elif hasattr(res.model, "lnalpha"):
            a = float(np.exp(res.model.lnalpha))
        else:
            # fallback: NaN
            a = np.nan
        return beta, se, pval, a
    except Exception:
        return np.nan, np.nan, 1.0, np.nan


def main():
    ap = argparse.ArgumentParser(
        prog="omicscanvas_nb_de.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            """Differential expression from raw counts using per-gene Negative Binomial regression (statsmodels).

This script performs size-factor normalization (DESeq median-of-ratios) and tests the treatment indicator coefficient.

Input must be RAW COUNTS (integers). Do NOT use FPKM/TPM as counts."""
        ),
        epilog=(
            """Examples:
  python omicscanvas_nb_de.py --counts counts.tsv --meta meta.tsv --control WT --treatment Mut --out de.tsv

counts.tsv format (genes x samples):
  gene_id	S1	S2	S3	S4
  GeneA	10	12	0	3

meta.tsv format:
  sample	condition
  S1	WT
  S2	WT
  S3	Mut
  S4	Mut
"""
        )
    )

    ap.add_argument("--counts", required=True, help="Raw counts matrix (genes x samples).")
    ap.add_argument("--meta", required=True, help="Sample metadata table with sample -> condition mapping.")
    ap.add_argument("--out", required=True, help="Output TSV for DE results.")

    ap.add_argument("--sample-col", default="sample", help="Column name for sample IDs in metadata.")
    ap.add_argument("--condition-col", default="condition", help="Column name for condition/group labels in metadata.")
    ap.add_argument("--control", required=True, help="Control condition name (reference group).")
    ap.add_argument("--treatment", required=True, help="Treatment condition name.")

    ap.add_argument(
        "--min-total-count", type=int, default=10,
        help="Prefilter: keep genes with total counts >= this threshold in the compared samples."
    )
    ap.add_argument(
        "--min-map-genes", type=int, default=10,
        help="Minimum number of genes required to compute size factors."
    )

    args = ap.parse_args()

    counts = read_counts(args.counts)
    meta = read_meta(args.meta, args.sample_col, args.condition_col)

    # Align samples
    common = counts.columns.intersection(meta.index)
    if len(common) < 4:
        eprint(f"[WARN] Only {len(common)} common samples between counts and meta.")
    if len(common) < 2:
        raise ValueError("No matching samples between counts columns and metadata sample IDs.")

    counts = counts[common]
    meta = meta.loc[common]

    # Keep only control + treatment
    cond = meta[args.condition_col].astype(str)
    keep_samples = cond.isin([args.control, args.treatment])
    meta = meta.loc[keep_samples]
    counts = counts[meta.index]  # match order

    # check replicates
    n_ctrl = int((meta[args.condition_col] == args.control).sum())
    n_trt = int((meta[args.condition_col] == args.treatment).sum())
    if n_ctrl < 2 or n_trt < 2:
        eprint(f"[WARN] Replicates are low (control={n_ctrl}, treatment={n_trt}). "
               "NB dispersion estimation may be unstable. Recommend >=3 per group.")

    # Prefilter
    gene_total = counts.sum(axis=1)
    counts = counts.loc[gene_total >= args.min_total_count]
    eprint(f"[INFO] Genes kept after prefilter: {counts.shape[0]}")

    if counts.shape[0] < args.min_map_genes:
        raise RuntimeError("Too few genes after filtering to compute size factors/model.")

    # Size factors (DESeq median ratio)
    sf = deseq_median_ratio_size_factors(counts)
    # offset is log(size_factor)
    offset = np.log(sf.values)

    # Design matrix: intercept + treatment indicator
    condition = meta[args.condition_col].astype(str).values
    x = np.column_stack([
        np.ones(len(condition), dtype=float),
        (condition == args.treatment).astype(float)
    ])

    # Fit per gene
    betas, ses, pvals, alphas, baseMeans = [], [], [], [], []

    # baseMean computed from normalized counts (counts / size_factor)
    norm = counts.div(sf, axis=1)
    base_mean = norm.mean(axis=1)

    y_mat = counts.to_numpy(dtype=float)  # genes x samples
    for i, gid in enumerate(counts.index):
        y = y_mat[i, :]
        beta, se, pval, alpha = fit_nb_for_gene(y, x, offset=offset)
        betas.append(beta)
        ses.append(se)
        pvals.append(pval if np.isfinite(pval) else 1.0)
        alphas.append(alpha)
        baseMeans.append(float(base_mean.loc[gid]))

        if (i + 1) % 2000 == 0:
            eprint(f"[INFO] Fitted {i+1} genes...")

    pvals = np.array(pvals, dtype=float)
    pvals = np.where(np.isnan(pvals), 1.0, pvals)
    padj = multipletests(pvals, method="fdr_bh")[1]

    # log2FC: beta is on log scale
    log2fc = np.array(betas, dtype=float) / np.log(2.0)

    res = pd.DataFrame({
        "gene_id": counts.index.astype(str),
        "baseMean": baseMeans,
        "log2FC": log2fc,
        "lfcSE": ses,
        "pvalue": pvals,
        "padj": padj,
        "alpha": alphas
    }).sort_values(["padj", "pvalue"])

    res.to_csv(args.out, sep="\t", index=False)
    eprint(f"[DONE] Results saved: {args.out}")


if __name__ == "__main__":
    main()
