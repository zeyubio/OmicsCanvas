# 08_compute_nb_de.py — Negative Binomial Differential Expression (NB-DE)

This script performs **pairwise differential expression** from **raw read counts** using **per-gene Negative Binomial regression** (statsmodels), with **DESeq-style median-of-ratios** size-factor normalization and **Benjamini–Hochberg FDR** correction.

It is intended as a lightweight DE tool when you already have a **counts matrix** and a **sample metadata table**.

---

## 1) What it does

Given:

- `counts.tsv`: genes × samples raw counts
- `meta.tsv`: sample → condition mapping
- one **control** condition and one **treatment** condition

The script:

1) Aligns samples between counts columns and metadata.
2) Keeps only samples in `{control, treatment}`.
3) Prefilters genes by total counts.
4) Computes **size factors** (DESeq median-of-ratios).
5) Fits per-gene NB regression:

- Design matrix: intercept + treatment indicator
- Offset: `log(size_factor)`  
- Tests the treatment coefficient

6) Outputs per-gene results and adjusts p-values by BH-FDR.

**Important:** input must be **RAW COUNTS** (integers). Do **NOT** use FPKM/TPM as counts.

---

## 2) Requirements

- Python ≥ 3.8
- `numpy`, `pandas`, `statsmodels`

Install:

```bash
pip install numpy pandas statsmodels
```

---

## 3) Inputs

### 3.1 Counts table (`--counts`) — required

Format: genes × samples.

- First column: `gene_id`
- Remaining columns: sample names (must match metadata sample IDs)

Example:

```tsv
gene_id  S1  S2  S3  S4
GeneA    10  12  0   3
GeneB    0   1   0   2
```

Rules:
- Values must be non-negative.
- Script expects integers; if non-integers are found, it warns and rounds.

### 3.2 Metadata table (`--meta`) — required

Tabular file with at least two columns:

- sample column (default `sample`)
- condition column (default `condition`)

Example:

```tsv
sample   condition
S1       WT
S2       WT
S3       Mut
S4       Mut
```

---

## 4) Output (`--out`)

A TSV with columns:

- `gene_id`
- `baseMean` : mean of normalized counts (counts / size_factor)
- `log2FC`   : estimated log2 fold-change (treatment vs control)
- `lfcSE`    : standard error of the treatment coefficient (on log2 scale)
- `pvalue`
- `padj`     : BH-FDR adjusted p-value
- `alpha`    : dispersion parameter (may be NaN if fitting fails)

Sorted by `padj`, then `pvalue`.

---

## 5) Usage examples

### 5.1 Basic DE: treatment vs control

```bash
python 08_compute_nb_de.py   --counts counts.tsv   --meta meta.tsv   --control WT   --treatment Mut   --out de.tsv
```

### 5.2 Custom metadata column names

```bash
python 08_compute_nb_de.py   --counts counts.tsv   --meta meta.tsv   --sample-col SampleID   --condition-col Group   --control WT   --treatment Mut   --out de.tsv
```

---

## 6) Arguments

Required:
- `--counts FILE` : raw counts matrix (genes × samples)
- `--meta FILE` : sample metadata
- `--control STR` : control condition name
- `--treatment STR` : treatment condition name
- `--out FILE` : output TSV

Metadata columns:
- `--sample-col STR` (default: `sample`)
- `--condition-col STR` (default: `condition`)

Filtering / robustness:
- `--min-total-count INT` (default: `10`)  
  Prefilter genes by total counts across the compared samples.
- `--min-map-genes INT` (default: `10`)  
  Minimum number of genes needed to compute size factors.

---

## 7) Notes & best practices

- Recommend ≥3 replicates per group for stable dispersion estimation. With <2 per group, NB inference is usually unreliable.
- Ensure sample names match exactly between the counts columns and metadata.
- If the output has many NaN dispersions / high p-values:
  - increase sequencing depth / use more replicates
  - relax `--min-total-count`
  - ensure you are using raw counts (not normalized)

---

## 8) Troubleshooting

### “Counts contain NA/non-numeric entries”
Your counts file has missing values or non-numeric tokens. Fix upstream.

### “No matching samples between counts columns and metadata sample IDs”
Check sample IDs and delimiters; they must match exactly.

### “Too few genes with positive counts to compute size factors”
Counts are too sparse after filtering; reduce `--min-total-count` or include more genes.

---

## 9) Help

```bash
python 08_compute_nb_de.py -h
```
