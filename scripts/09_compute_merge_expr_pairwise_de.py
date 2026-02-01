#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
omicscanvas_merge_expr_and_pairwise_de.py

Purpose
-------
This is an end-to-end helper for expression matrices and pairwise DE based on a directory of per-sample quantification tables.

It:
1) Scans an input directory and collects sample tables (via --pattern / --suffix / --auto-discover).
2) Infers or reads sample -> condition mapping (replicates) from filenames (regex) or a metadata file.
3) Writes 4 base outputs:
   - merged_counts.tsv
   - merged_FPKM.tsv
   - merged_TPM.tsv
   - sample_meta.tsv
4) Runs pairwise comparisons:
   - By default: only A_vs_B (not both directions). If --control is provided: all conditions vs control.
   - For each pair: keep gene_id, meanA, meanB, log2(meanA/meanB), pvalue, padj, alpha.
5) Produces two final summary files:
   - final_FPKM_summary.tsv
   - final_TPM_summary.tsv

Statistics
----------
- pvalue/padj/alpha are estimated from raw counts using per-gene Negative Binomial regression (statsmodels).
- log2FC is computed from mean(FPKM) or mean(TPM) with a pseudocount.
- Low-expression filtering (min_total_count) affects NB pvalue/padj/alpha only.

Dependencies
------------
  pip install pandas numpy statsmodels
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

from statsmodels.discrete.discrete_model import NegativeBinomial
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning, HessianInversionWarning
from statsmodels.stats.multitest import multipletests


# --------------------------
# basic helpers
# --------------------------
def eprint(msg: str) -> None:
    print(msg, file=sys.stderr)


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def safe_token(x: str) -> str:
    return "".join(c if (c.isalnum() or c in "._-+") else "_" for c in x)


def sniff_has_required_cols(path: Path, required: List[str]) -> bool:
    try:
        df = pd.read_csv(path, sep=None, engine="python", nrows=5)
    except Exception:
        return False
    cols = set(map(str, df.columns))
    return all(c in cols for c in required)


def discover_files(indir: Path,
                   *,
                   suffix: Optional[str],
                   pattern: Optional[str],
                   auto_discover: bool,
                   exts: List[str],
                   required_cols: List[str]) -> List[Path]:
    # 1) --pattern has highest priority
    if pattern:
        files = sorted([p for p in indir.glob(pattern) if p.is_file()])
        if files:
            return files

    # 2) suffix
    if suffix:
        files = sorted([p for p in indir.iterdir() if p.is_file() and p.name.endswith(suffix)])
        if files:
            return files

    # 3) auto discover
    if not auto_discover:
        raise FileNotFoundError(
            "No input files found. Use --pattern or --suffix, or enable --auto-discover."
        )

    cand = []
    for p in indir.iterdir():
        if not p.is_file():
            continue
        low = p.name.lower()
        if any(low.endswith(e) for e in exts):
            cand.append(p)
    cand = sorted(cand)
    files = [p for p in cand if sniff_has_required_cols(p, required_cols)]
    if not files:
        raise FileNotFoundError(
            f"Auto-discovery failed: no text file in {indir} contains required columns {required_cols}.\n"
            f"Use --pattern or --suffix to explicitly specify input files."
        )
    return files


def sample_base_from_filename(p: Path, suffix: Optional[str]) -> str:
    name = p.name
    if suffix and name.endswith(suffix):
        return name[:-len(suffix)]
    return p.stem


def infer_condition_by_regex(sample_base: str, rep_regex: str) -> str:
    import re
    m = re.match(rep_regex, sample_base)
    if m:
        return str(m.group(1))
    return sample_base


def infer_condition_auto(sample_base: str) -> str:
    import re
    patterns = [
        r"^(.+?)(?:[_-]\d+)$",                  # _1 / -2
        r"^(.+?)(?:[_-]rep\d+)$",               # _rep1
        r"^(.+?)(?:[_-]r\d+)$",                 # _r1
        r"^(.+?)(?:[._-]rep\d+)$",              # .rep1
        r"^(.+?)(?:[._-]R\d+)$",                # .R1
        r"^(.+?)(?:[_-]R\d+)$",                 # _R2
    ]
    for pat in patterns:
        m = re.match(pat, sample_base)
        if m:
            return str(m.group(1))
    return sample_base


def read_one_expr_table(path: Path,
                        id_col: str,
                        counts_col: str,
                        fpkm_col: str,
                        tpm_col: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    for col in (id_col, counts_col, fpkm_col, tpm_col):
        if col not in df.columns:
            raise ValueError(f"Missing column '{col}' in {path.name}. Columns: {list(df.columns)[:50]}")
    df = df[[id_col, counts_col, fpkm_col, tpm_col]].copy()
    df[id_col] = df[id_col].astype(str)
    for col in (counts_col, fpkm_col, tpm_col):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if df[[counts_col, fpkm_col, tpm_col]].isna().any().any():
        bad = int(df[[counts_col, fpkm_col, tpm_col]].isna().sum().sum())
        raise ValueError(f"NA/non-numeric entries detected in {path.name}: {bad}")
    if (df[counts_col] < 0).any():
        raise ValueError(f"Negative counts detected in {path.name}")
    return df


# --------------------------
# NB DE core (counts)
# --------------------------
def deseq_median_ratio_size_factors(counts: pd.DataFrame, min_map_genes: int = 10) -> pd.Series:
    X = counts.to_numpy(dtype=float)  # genes x samples
    with np.errstate(divide="ignore", invalid="ignore"):
        logX = np.where(X > 0, np.log(X), np.nan)
    gm = np.exp(np.nanmean(logX, axis=1))
    valid_gene = np.isfinite(gm) & (gm > 0)

    if int(valid_gene.sum()) < int(min_map_genes):
        raise RuntimeError("Too few genes with positive counts to compute size factors.")

    ratios = X[valid_gene, :] / gm[valid_gene, None]
    sf = []
    for j in range(ratios.shape[1]):
        r = ratios[:, j]
        r = r[np.isfinite(r) & (r > 0)]
        sf.append(np.median(r) if len(r) else np.nan)
    sf = np.array(sf, dtype=float)

    if (not np.all(np.isfinite(sf))) or np.any(sf <= 0):
        raise RuntimeError("Failed to compute valid size factors. Check counts matrix.")
    sf = sf / np.exp(np.mean(np.log(sf)))
    return pd.Series(sf, index=counts.columns, name="size_factor")


def fit_nb_for_gene(y: np.ndarray, x: np.ndarray, offset: np.ndarray) -> Tuple[float, float, float]:
    """
    Return: pvalue, alpha, ok_flag (0/1 as float)

    Note: statsmodels.discrete.NegativeBinomial params = [beta..., alpha] (alpha is already positive; not lnalpha).
    """
    try:
        model = NegativeBinomial(y, x, offset=offset)

        # Suppress per-gene fitting warnings (we use ok_flag to mark failed/unreliable fits)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = model.fit(disp=False, maxiter=200)

        converged = True
        if hasattr(res, "mle_retvals") and isinstance(res.mle_retvals, dict):
            converged = bool(res.mle_retvals.get("converged", True))

        # If ConvergenceWarning occurs, treat as unreliable convergence
        for ww in w:
            if issubclass(ww.category, ConvergenceWarning):
                converged = False
                break

        # pvalue: if Hessian inversion fails, pvalues may be NaN; set to 1.0
        pval = float(res.pvalues[1]) if len(res.pvalues) > 1 else 1.0
        if (not np.isfinite(pval)) or (not converged):
            pval = 1.0

        # alpha: typically the last parameter
        alpha = np.nan
        if hasattr(res, "params") and len(res.params) == (x.shape[1] + 1):
            alpha = float(res.params[-1])
        elif hasattr(res, "params") and len(res.params) >= 1:
            # Fallback: use the last parameter as alpha (rare edge cases across versions)
            alpha = float(res.params[-1])

        ok = 1.0 if (converged and np.isfinite(alpha)) else 0.0
        return pval, alpha, ok

    except Exception:
        return 1.0, np.nan, 0.0


def nb_pvals_for_pair(counts: pd.DataFrame,
                      cond: pd.Series,
                      control: str,
                      treatment: str,
                      min_total_count: int,
                      min_map_genes: int) -> pd.DataFrame:
    """
    Return a DataFrame indexed by gene_id:
      pvalue, padj, alpha
    Note: only fit genes passing filters; others remain NA (can be joined later).
    Filters: total counts >= min_total_count, and mean counts in each group must be > 0.
    """
    keep = cond.isin([control, treatment])
    cond2 = cond.loc[keep].astype(str)
    counts2 = counts.loc[:, cond2.index]

    n_ctrl = int((cond2 == control).sum())
    n_trt = int((cond2 == treatment).sum())
    if n_ctrl < 2 or n_trt < 2:
        eprint(f"[WARN] Replicates low for NB: {control}={n_ctrl}, {treatment}={n_trt} (recommend >=3).")

    # Prefilter: genes with mean counts == 0 in either group (reduces alpha/Hessian failures)
    gene_total = counts2.sum(axis=1)
    mean_ctrl = counts2.loc[:, (cond2 == control)].mean(axis=1)
    mean_trt  = counts2.loc[:, (cond2 == treatment)].mean(axis=1)

    use = (gene_total >= int(min_total_count)) & (mean_ctrl > 0) & (mean_trt > 0)

    # Record filter rate (helps diagnose excessive NB failures)
    dropped_zero_mean = int(((mean_ctrl <= 0) | (mean_trt <= 0)).sum())
    if dropped_zero_mean > 0:
        eprint(f"[INFO] Drop {dropped_zero_mean} genes with zero mean counts in either group: {treatment} vs {control}")

    counts_f = counts2.loc[use].copy()

    if counts_f.shape[0] < int(min_map_genes):
        eprint("[WARN] Too few genes after filtering; pvalue/padj will be NA.")
        return pd.DataFrame({"pvalue": pd.Series(dtype=float),
                             "padj": pd.Series(dtype=float),
                             "alpha": pd.Series(dtype=float)})

    if not np.all(np.mod(counts_f.values, 1) == 0):
        eprint("[WARN] Counts contain non-integers; rounding to nearest int.")
    counts_f = counts_f.round().astype(int)

    sf = deseq_median_ratio_size_factors(counts_f, min_map_genes=min_map_genes)
    offset = np.log(sf.values)

    x = np.column_stack([
        np.ones(len(cond2), dtype=float),
        (cond2.values == treatment).astype(float)
    ])

    y_mat = counts_f.to_numpy(dtype=float)
    gids = counts_f.index.astype(str).tolist()

    pvals = np.empty(len(gids), dtype=float)
    alphas = np.empty(len(gids), dtype=float)

    for i, gid in enumerate(gids):
        y = y_mat[i, :]
        pval, alpha, ok = fit_nb_for_gene(y, x, offset=offset)
        pvals[i] = pval
        alphas[i] = alpha
        if (i + 1) % 5000 == 0:
            eprint(f"[INFO] NB fitted {i+1} genes for {treatment} vs {control} ...")

    padj = multipletests(pvals, method="fdr_bh")[1]

    out = pd.DataFrame({"pvalue": pvals, "padj": padj, "alpha": alphas}, index=gids)
    out.index.name = "gene_id"
    return out


# --------------------------
# assembling outputs
# --------------------------
def mean_by_condition(mat: pd.DataFrame, cond: pd.Series, conditions: List[str]) -> pd.DataFrame:
    out = {}
    for c in conditions:
        cols = cond[cond == c].index.tolist()
        out[c] = mat[cols].mean(axis=1) if cols else np.nan
    df = pd.DataFrame(out)
    df.index = df.index.astype(str)
    return df


def build_pair_minimal(metric_means: pd.DataFrame,
                       nb_stats: pd.DataFrame,
                       a1: str,
                       a2: str,
                       pseudocount: float) -> pd.DataFrame:
    """
    metric_means: genes x conditions (mean matrix, e.g., mean_fpkm / mean_tpm)
    nb_stats: index=gene_id, columns pvalue/padj/alpha (only for fitted genes)
    Outputs: gene_id, a1_mean, a2_mean, log2(a1/a2), pvalue, padj, alpha
    """
    df = pd.DataFrame({
        "gene_id": metric_means.index.astype(str),
        "a1_mean": metric_means[a1].values,
        "a2_mean": metric_means[a2].values,
    })
    df["log2FC"] = np.log2((df["a1_mean"].astype(float) + pseudocount) / (df["a2_mean"].astype(float) + pseudocount))

    if nb_stats is not None and len(nb_stats) > 0:
        tmp = nb_stats.reset_index()
        df = df.merge(tmp, how="left", on="gene_id")
    else:
        df["pvalue"] = np.nan
        df["padj"] = np.nan
        df["alpha"] = np.nan

    # Keep column order
    df = df[["gene_id", "a1_mean", "a2_mean", "log2FC", "pvalue", "padj", "alpha"]]
    return df


def add_pair_cols_to_summary(summary: pd.DataFrame,
                             pair_df: pd.DataFrame,
                             *,
                             pair_name: str) -> pd.DataFrame:
    """
    Add to summary:
      log2FC_<pair>, pvalue_<pair>, padj_<pair>, alpha_<pair>
    """
    pair_df = pair_df.set_index("gene_id")
    summary[f"log2FC_{pair_name}"] = pair_df["log2FC"]
    summary[f"pvalue_{pair_name}"] = pair_df["pvalue"]
    summary[f"padj_{pair_name}"] = pair_df["padj"]
    summary[f"alpha_{pair_name}"] = pair_df["alpha"]
    return summary


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Merge expression tables from a folder, infer replicates/conditions, and run pairwise NB DE (minimal outputs)."
    )
    ap.add_argument("--indir", required=True, help="Input folder containing expression tables.")
    ap.add_argument("--outdir", required=True, help="Output folder.")

    # file discovery
    ap.add_argument("--suffix", default=None, help="Match files by suffix (e.g. _reads_count.txt).")
    ap.add_argument("--pattern", default=None, help="Glob pattern under indir (e.g. '*_reads_count.txt'). Overrides --suffix.")
    ap.add_argument("--auto-discover", action="store_true", help="If pattern/suffix finds nothing, scan common text files and auto-detect by required columns.")
    ap.add_argument("--auto-exts", default=".txt,.tsv,.csv", help="Extensions for auto-discover (comma-separated).")
    ap.add_argument("--require-cols", default="ID,FPKM,TPM,counts", help="Columns required for auto-discover (comma-separated).")

    # columns
    ap.add_argument("--id-col", default="ID")
    ap.add_argument("--counts-col", default="counts")
    ap.add_argument("--fpkm-col", default="FPKM")
    ap.add_argument("--tpm-col", default="TPM")

    # replicate parsing
    ap.add_argument("--rep-auto", action="store_true", help="Auto-infer condition by stripping common replicate tokens.")
    ap.add_argument("--rep-regex", default=r"(.+?)(?:_rep\d+|_\d+)$",
                    help="Custom regex stripping replicate suffix; must capture condition in group(1). Ignored if --rep-auto is set.")

    # comparisons control
    ap.add_argument("--wt", default=None,
                    help="Optional WT/control condition name(s), comma-separated. If set, compare all non-WT vs WT.")
    ap.add_argument("--wt-merged-name", default="WT",
                    help="If multiple WT names are provided, merge them into this control name.")

    # DE parameters
    ap.add_argument("--min-total-count", type=int, default=10,
                    help="NB prefilter genes: keep genes with total counts >= threshold in compared samples.")
    ap.add_argument("--min-map-genes", type=int, default=10,
                    help="Minimum genes for size factor estimation.")
    ap.add_argument("--pseudocount", type=float, default=1e-6,
                    help="Pseudocount for log2(mean ratio) for FPKM/TPM.")

    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    ensure_dir(outdir)
    ensure_dir(outdir / "matrices")
    ensure_dir(outdir / "pairs_minimal")
    ensure_dir(outdir / "final")

    required_cols = [x.strip() for x in args.require_cols.split(",") if x.strip()]
    exts = [x.strip() for x in args.auto_exts.split(",") if x.strip()]
    if not exts:
        exts = [".txt", ".tsv", ".csv"]

    files = discover_files(
        indir,
        suffix=args.suffix,
        pattern=args.pattern,
        auto_discover=args.auto_discover,
        exts=exts,
        required_cols=required_cols
    )
    eprint(f"[INFO] Found {len(files)} input files.")

    counts_dict: Dict[str, pd.Series] = {}
    fpkm_dict: Dict[str, pd.Series] = {}
    tpm_dict: Dict[str, pd.Series] = {}
    sample2cond: Dict[str, str] = {}

    for p in files:
        sample_base = sample_base_from_filename(p, args.suffix)
        if args.rep_auto:
            cond = infer_condition_auto(sample_base)
        else:
            cond = infer_condition_by_regex(sample_base, args.rep_regex)
        if sample_base in sample2cond:
            raise ValueError(f"Duplicate sample name after parsing: {sample_base}")
        sample2cond[sample_base] = cond

        df = read_one_expr_table(p, args.id_col, args.counts_col, args.fpkm_col, args.tpm_col).set_index(args.id_col)

        counts_dict[sample_base] = df[args.counts_col].astype(float)
        fpkm_dict[sample_base] = df[args.fpkm_col].astype(float)
        tpm_dict[sample_base] = df[args.tpm_col].astype(float)

    counts = pd.DataFrame(counts_dict).fillna(0.0)
    fpkm = pd.DataFrame(fpkm_dict).fillna(0.0)
    tpm = pd.DataFrame(tpm_dict).fillna(0.0)

    # counts -> int
    if not np.all(np.mod(counts.values, 1) == 0):
        eprint("[WARN] Some counts are non-integers; rounding.")
    counts = counts.round().astype(int)

    # write required matrices
    counts.to_csv(outdir / "matrices" / "merged_counts.tsv", sep="\t")
    fpkm.to_csv(outdir / "matrices" / "merged_FPKM.tsv", sep="\t")
    tpm.to_csv(outdir / "matrices" / "merged_TPM.tsv", sep="\t")

    meta = pd.DataFrame({"sample": list(sample2cond.keys()), "condition": list(sample2cond.values())}).set_index("sample")
    meta.to_csv(outdir / "matrices" / "sample_meta.tsv", sep="\t")

    cond_series = meta["condition"].copy()

    # merge WT names if provided
    if args.wt:
        wt_names = [x.strip() for x in args.wt.split(",") if x.strip()]
        if not wt_names:
            raise ValueError("--wt provided but empty after parsing.")
        cond_series = cond_series.where(~cond_series.isin(wt_names), other=args.wt_merged_name)

    conditions = sorted(cond_series.unique().tolist())
    eprint(f"[INFO] Conditions: {conditions}")

    # compute means
    mean_counts = mean_by_condition(counts, cond_series, conditions)
    mean_fpkm = mean_by_condition(fpkm, cond_series, conditions)
    mean_tpm = mean_by_condition(tpm, cond_series, conditions)

    # define comparisons (ONLY one direction)
    comparisons: List[Tuple[str, str]] = []  # (a1, a2) meaning log2(a1/a2) and NB treatment=a1, control=a2
    if args.wt:
        ctrl = args.wt_merged_name
        if ctrl not in conditions:
            raise ValueError(f"WT merged name '{ctrl}' not found in inferred conditions: {conditions}")
        for c in conditions:
            if c == ctrl:
                continue
            comparisons.append((c, ctrl))
    else:
        # unordered pairs once (A_vs_B where A is later, B is earlier, deterministic)
        for i in range(len(conditions)):
            for j in range(i + 1, len(conditions)):
                a2 = conditions[i]
                a1 = conditions[j]
                comparisons.append((a1, a2))

    eprint(f"[INFO] Total comparisons: {len(comparisons)}")

    # initialize final summary tables:
    # ID + mean columns
    fpkm_summary = mean_fpkm.copy()
    fpkm_summary.insert(0, "ID", fpkm_summary.index.astype(str))
    tpm_summary = mean_tpm.copy()
    tpm_summary.insert(0, "ID", tpm_summary.index.astype(str))

    # run pairs and add columns
    for a1, a2 in comparisons:
        pair_name = f"{safe_token(a1)}_vs_{safe_token(a2)}"
        eprint(f"[INFO] Pair: {a1} vs {a2}")

        # NB stats from counts
        nb_stats = nb_pvals_for_pair(
            counts=counts,
            cond=cond_series,
            control=a2,
            treatment=a1,
            min_total_count=args.min_total_count,
            min_map_genes=args.min_map_genes
        )

        # minimal pair outputs for FPKM and TPM (same pvalue/padj/alpha; different means/log2FC)
        pair_fpkm = build_pair_minimal(mean_fpkm, nb_stats, a1, a2, args.pseudocount)
        pair_tpm  = build_pair_minimal(mean_tpm,  nb_stats, a1, a2, args.pseudocount)

        # write pair files
        pair_fpkm.to_csv(outdir / "pairs_minimal" / f"{pair_name}.FPKM.tsv", sep="\t", index=False)
        pair_tpm.to_csv(outdir / "pairs_minimal" / f"{pair_name}.TPM.tsv", sep="\t", index=False)

        # add columns into final summary
        fpkm_summary = add_pair_cols_to_summary(
            fpkm_summary.set_index("ID"),
            pair_fpkm,
            pair_name=pair_name
        ).reset_index().rename(columns={"index": "ID"} if "index" in fpkm_summary.columns else {})

        tpm_summary = add_pair_cols_to_summary(
            tpm_summary.set_index("ID"),
            pair_tpm,
            pair_name=pair_name
        ).reset_index().rename(columns={"index": "ID"} if "index" in tpm_summary.columns else {})

    # final outputs
    fpkm_summary.to_csv(outdir / "final" / "final_FPKM_summary.tsv", sep="\t", index=False)
    tpm_summary.to_csv(outdir / "final" / "final_TPM_summary.tsv", sep="\t", index=False)

    eprint("[DONE] Outputs written:")
    eprint(f"  {outdir}/matrices/merged_counts.tsv")
    eprint(f"  {outdir}/matrices/merged_FPKM.tsv")
    eprint(f"  {outdir}/matrices/merged_TPM.tsv")
    eprint(f"  {outdir}/matrices/sample_meta.tsv")
    eprint(f"  {outdir}/pairs_minimal/*.(FPKM|TPM).tsv")
    eprint(f"  {outdir}/final/final_FPKM_summary.tsv")
    eprint(f"  {outdir}/final/final_TPM_summary.tsv")


if __name__ == "__main__":
    main()
