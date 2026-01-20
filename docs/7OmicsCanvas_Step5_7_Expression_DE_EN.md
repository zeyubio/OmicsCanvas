# OmicsCanvas | Expression Quantification + Differential Expression (Step 5–7)

This document describes three command-line utilities (renamed to remove numeric prefixes):

- **Step 5**: `omicscanvas_bam_to_fpkm.py` — quantify gene-level expression from a BAM (FPKM/TPM + raw counts).
- **Step 6**: `omicscanvas_nb_de.py` — differential expression from a merged **raw counts** matrix via Negative Binomial regression.
- **Step 7**: `omicscanvas_merge_expr_and_pairwise_de.py` — merge many per-sample expression files from a directory and run pairwise DE (and summary tables).

---

## Step 5: omicscanvas_bam_to_fpkm.py

### What it does
Given:
- a gene BED file (genomic intervals),
- a gene length table (bp), and
- a coordinate-sorted, indexed BAM,

the script outputs a TSV with per-gene:
- `counts` (strict read/fragment counts in the interval),
- `FPKM`,
- `TPM`,
- `length_bp`.

### Important definitions
- **Single-end (SE)**: counts reads passing filters in each interval.
- **Paired-end (PE)**: counts **fragments**, deduplicated by `query_name` **within each interval**.
- **Library size (N)**: computed using the *same* strict filters as used for interval counting.

> If gene intervals overlap, the same read/fragment may be counted for multiple genes. This is expected for interval-based counting.

### Input formats
**BED** (minimum 4 columns, tab-separated):
```
chrom   start   end     gene_id
Chr01   100     500     GeneA
```
- By default, BED is **0-based, end-exclusive**.
- If your BED is **1-based, inclusive**, use `--bed-one-based`.

**Length table** (2 columns, tab-separated, no header):
```
GeneA   1200
GeneB   980
```

### Command examples
**Single-end**:
```bash
python omicscanvas_bam_to_fpkm.py   --bed genes.bed   --length gene_len.tsv   --bam sample.bam   --out sample.expr.tsv   --mode se   --min-mapq 10
```

**Paired-end (strict proper pairs)**:
```bash
python omicscanvas_bam_to_fpkm.py   --bed genes.bed   --length gene_len.tsv   --bam sample.bam   --out sample.expr.tsv   --mode pe   --min-mapq 10   --threads 8
```

### Parameters
- `-b/--bed` (str, required): gene intervals BED.
- `-l/--length` (str, required): gene length table (bp).
- `-m/--bam` (str, required): coordinate-sorted + indexed BAM (`.bai`).
- `-o/--out` (str, required): output TSV.
- `--mode` (`auto|se|pe`, default `auto`): read layout.
- `--min-mapq` (int, default `0`): minimum MAPQ.
- `--include-duplicates` (flag): include duplicate reads.
- `--include-qcfail` (flag): include QC-failed reads.
- `--require-proper-pair` (flag, PE): require proper pairs (recommended).
- `--allow-improper-pair` (flag, PE): allow improper pairs (disables proper-pair requirement).
- `--bed-one-based` (flag): convert 1-based inclusive BED to 0-based half-open.
- `--threads` (int, default `1`): pysam BGZF decompression threads.
- `--progress` (int, default `1000`): print progress every N intervals (`0` disables).

---

## Step 6: omicscanvas_nb_de.py

### What it does
Per-gene **Negative Binomial (NB2) regression** on **raw counts**:
- size-factor normalization: DESeq median-of-ratios,
- design: `y ~ intercept + I(condition==treatment)`,
- outputs: `log2FC`, `pvalue`, `padj` (BH-FDR), plus `baseMean`, `lfcSE`, `alpha`.

### Inputs
1) **Counts table** (genes × samples; first column gene_id):
```
gene_id S1 S2 S3 S4
GeneA   10 12 0  3
```
2) **Metadata** mapping sample → condition:
```
sample condition
S1     WT
S2     WT
S3     Mut
S4     Mut
```

### Command example
```bash
python omicscanvas_nb_de.py   --counts counts.tsv   --meta meta.tsv   --sample-col sample   --condition-col condition   --control WT   --treatment Mut   --out de.tsv
```

### Parameters
- `--counts` (str, required): raw counts matrix.
- `--meta` (str, required): metadata table.
- `--out` (str, required): output TSV.
- `--sample-col` (str, default `sample`): column name for sample IDs.
- `--condition-col` (str, default `condition`): column name for group labels.
- `--control` (str, required): control group label.
- `--treatment` (str, required): treatment group label.
- `--min-total-count` (int, default `10`): prefilter genes with total counts < threshold.
- `--min-map-genes` (int, default `10`): minimum genes needed to compute size factors.

### Notes
- NB dispersion estimation is unstable with very low replicates; **>=3 per group is recommended**.

---

## Step 7: omicscanvas_merge_expr_and_pairwise_de.py

### What it does
From a directory of per-sample expression files (counts/FPKM/TPM in one file), this script:
1) discovers input files (by `--pattern`, `--suffix`, or `--auto-discover`),
2) infers condition/replicates from file names (or you can override),
3) merges to:
   - `merged_counts.tsv`
   - `merged_FPKM.tsv`
   - `merged_TPM.tsv`
   - `sample_meta.tsv`
4) runs pairwise DE (NB on counts) and writes per-pair results,
5) produces summary tables:
   - `final_FPKM_summary.tsv`
   - `final_TPM_summary.tsv`

### Typical file naming
If your input files look like:
- `DK_col_rep1.expr.tsv`
- `DK_col_rep2.expr.tsv`
- `BL_col_rep1.expr.tsv`

then you can use a replicate-regex to collapse replicates:
- condition = `DK_col` / `BL_col`
- replicate = `rep1` / `rep2`

### Command example
```bash
python omicscanvas_merge_expr_and_pairwise_de.py   --indir expr_dir   --suffix .expr.tsv   --rep-regex "(.*)_rep\d+$"   --metric FPKM   --outdir out   --wt DK_col
```

### Parameters (high-level)
The script is feature-rich; run `--help` to see all options. The most important knobs:

- Input discovery:
  - `--indir` (required): input directory
  - `--pattern` (glob): strongest priority
  - `--suffix` (string): match file name suffix
  - `--auto-discover` (flag): attempt to detect suitable tables

- Column names (if your tables differ):
  - `--id-col` (default `ID`)
  - `--count-col` (default `counts`)
  - `--fpkm-col` (default `FPKM`)
  - `--tpm-col` (default `TPM`)

- Replicate parsing:
  - `--rep-regex`: regex to map sample file base name → condition

- DE options:
  - `--wt`: if set, compare every condition vs WT only
  - `--min-total-count`: prefilter for NB testing
  - `--pseudocount`: used in log2 ratio for summary tables

---

## Recommended workflow

### Option A (simple two-group DE)
1) Generate raw counts matrix (your own pipeline or Step 7 merged_counts.tsv)
2) Run Step 6 (`omicscanvas_nb_de.py`) for one comparison.

### Option B (full multi-condition pipeline)
1) For each sample BAM, run Step 5 to get `*.expr.tsv` (FPKM/TPM + counts).
2) Put all `*.expr.tsv` into one directory.
3) Run Step 7 to merge + run pairwise DE + export summary tables.

---

## Dependencies
- Step 5: `pysam`, `pandas`
- Step 6/7: `pandas`, `numpy`, `statsmodels`

