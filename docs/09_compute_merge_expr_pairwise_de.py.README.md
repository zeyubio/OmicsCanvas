# 09_compute_merge_expr_pairwise_de.py — Merge Expression Tables + Pairwise NB-DE

This script is an end-to-end helper that:

1) Scans a directory of **per-sample expression tables**
2) Merges them into **counts / FPKM / TPM matrices**
3) Infers **sample → condition** (replicates) from filenames (auto or regex)
4) Runs **pairwise differential expression** (NB regression on counts)
5) Outputs:
   - minimal pairwise tables (means + log2FC + pvalue/padj/alpha)
   - final summary tables with pairwise columns appended

NB statistics (`pvalue/padj/alpha`) are computed from **raw counts** only.
`log2FC` is computed from mean(FPKM) or mean(TPM) using a **pseudocount**.

---

## 1) Requirements

- Python ≥ 3.8
- `numpy`, `pandas`, `statsmodels`

Install:

```bash
pip install numpy pandas statsmodels
```

---

## 2) Input files (per-sample expression tables)

Each file must contain the following columns (defaults shown; configurable):

- `ID` (gene_id)
- `counts` (raw counts)
- `FPKM`
- `TPM`

You can change column names via:
- `--id-col`
- `--counts-col`
- `--fpkm-col`
- `--tpm-col`

---

## 3) File discovery

The script discovers per-sample files in `--indir` by priority:

1) `--pattern` (glob) if provided and matches anything
2) `--suffix` if provided and matches anything
3) If still nothing and `--auto-discover` is set:
   - scan text-like files in `--auto-exts`
   - keep files that contain all `--require-cols`

Recommended:
- If your filenames follow a stable pattern, use `--suffix` or `--pattern`.

---

## 4) Sample names and condition (replicate) inference

For each input file:
- sample name is derived from filename (stem), optionally removing `--suffix`
- condition is inferred by either:

A) `--rep-auto`  
   strips common replicate tokens like `_1`, `_rep1`, `_r1`, `.R1`, etc.

B) `--rep-regex` (default)  
   must capture the condition in **group(1)** (e.g. `(.+?)(?:_rep\d+|_\d+)$`)

The script writes the inferred mapping to:
- `<outdir>/matrices/sample_meta.tsv`

---

## 5) Comparisons

Two modes:

### 5.1 Control mode (`--wt`)
If `--wt` is provided (comma-separated condition names), the script merges them into `--wt-merged-name` (default `WT`), then runs:

- each non-WT condition vs WT

### 5.2 All-pairs mode (default)
If `--wt` is not provided, the script runs each unordered pair **once** (deterministic direction):
- later condition vs earlier condition (alphabetical order of condition names)

---

## 6) Outputs (folder structure)

Given `--outdir OUT`:

### 6.1 Matrices
- `OUT/matrices/merged_counts.tsv`
- `OUT/matrices/merged_FPKM.tsv`
- `OUT/matrices/merged_TPM.tsv`
- `OUT/matrices/sample_meta.tsv`

### 6.2 Pairwise minimal outputs
For each pair `<A>_vs_<B>`:
- `OUT/pairs_minimal/<A>_vs_<B>.FPKM.tsv`
- `OUT/pairs_minimal/<A>_vs_<B>.TPM.tsv`

Columns:
- `gene_id`
- `a1_mean` (mean of A)
- `a2_mean` (mean of B)
- `log2FC = log2((a1_mean+p)/(a2_mean+p))`
- `pvalue`, `padj`, `alpha` (from NB on counts)

### 6.3 Final summary tables
- `OUT/final/final_FPKM_summary.tsv`
- `OUT/final/final_TPM_summary.tsv`

They contain:
- per-condition mean columns
- plus appended columns for each pair:
  - `log2FC_<pair>`
  - `pvalue_<pair>`
  - `padj_<pair>`
  - `alpha_<pair>`

---

## 7) Usage examples

### 7.1 Simple: suffix-based discovery, all pairs
```bash
python 09_compute_merge_expr_pairwise_de.py   --indir quant_tables   --outdir out_pairwise   --suffix _reads_count.txt   --rep-auto
```

### 7.2 Control mode: compare all vs WT
```bash
python 09_compute_merge_expr_pairwise_de.py   --indir quant_tables   --outdir out_vs_wt   --suffix _reads_count.txt   --rep-auto   --wt WT   --wt-merged-name WT
```

### 7.3 Custom replicate parsing regex
```bash
python 09_compute_merge_expr_pairwise_de.py   --indir quant_tables   --outdir out_pairwise   --pattern "*_reads_count.txt"   --rep-regex "(.+?)_R\d+$"
```

### 7.4 Auto-discover when you don't know filenames
```bash
python 09_compute_merge_expr_pairwise_de.py   --indir quant_tables   --outdir out_pairwise   --auto-discover   --require-cols "ID,FPKM,TPM,counts"
```

---

## 8) Arguments

Required:
- `--indir DIR` : folder containing per-sample expression tables
- `--outdir DIR` : output folder

File discovery:
- `--suffix STR` : match files by suffix
- `--pattern STR` : glob pattern under indir (overrides --suffix)
- `--auto-discover` : enable auto discovery if suffix/pattern find nothing
- `--auto-exts STR` : extensions list for auto-discover (default: `.txt,.tsv,.csv`)
- `--require-cols STR` : required columns for auto-discover (default: `ID,FPKM,TPM,counts`)

Columns:
- `--id-col STR` (default: `ID`)
- `--counts-col STR` (default: `counts`)
- `--fpkm-col STR` (default: `FPKM`)
- `--tpm-col STR` (default: `TPM`)

Replicate inference:
- `--rep-auto` : auto infer condition name
- `--rep-regex STR` : regex capturing condition in group(1)

Control:
- `--wt STR` : comma-separated control condition name(s)
- `--wt-merged-name STR` (default: `WT`) : merged control name

NB-DE parameters:
- `--min-total-count INT` (default: `10`)
- `--min-map-genes INT` (default: `10`)
- `--pseudocount FLOAT` (default: `1e-6`) : for log2FC on mean(FPKM)/mean(TPM)

---

## 9) Notes & best practices

- NB-DE requires replicates. With <2 per condition, p-values can be unstable.
- If you see many NA p-values, try increasing `--min-total-count` or ensure each group has enough nonzero counts.
- Ensure your `counts` column is raw counts (not normalized).

---

## 10) Help

```bash
python 09_compute_merge_expr_pairwise_de.py -h
```
