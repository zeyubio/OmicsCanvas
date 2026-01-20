# OmicsCanvas


OmicsCanvas is a step-based pipeline for gene-centric multi-omics visualization and analysis (ChIP-seq/ATAC-seq/WGBS methylation, plus RNA expression).

It is organized into 16 numbered steps. **Steps 1–4** are preparation, **5–9** are computation, and **10–16** are plotting/analysis.


## Installation


### Option A: Portable (unzip & run)


1. Install Python 3.8+ and dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Run any step script in `scripts/`, or use the runner:
   ```bash
   python omicscanvas.py --list
   python omicscanvas.py 10 -- --help
   ```


### Option B: Pip install (from repo)


If you use the pip package version, you can do:
```bash
pip install .
omicscanvas --list
omicscanvas run 10 -- --help
```


## Pipeline overview


| Step | Phase | Purpose | Script |

|---:|---|---|---|

| 01 | prepare | GFF3/GTF -> BED + gene length table | `scripts/01_prepare_gff_to_bed_genes_length.py` |

| 02 | prepare | Split CX report by context (CG/CHG/CHH) | `scripts/02_prepare_cx_context_split.py` |

| 03 | prepare | Merge CX replicates | `scripts/03_prepare_cx_replicate_merge.py` |

| 04 | prepare | Extract per-gene methylation table from CX | `scripts/04_prepare_extract_gene_methylation.py` |

| 05 | compute | BAM -> gene matrix (TSS/gene/TES) | `scripts/05_compute_bam_to_gene_matrices.py` |

| 06 | compute | Build CX_gene matrices/profiles from CX/genes | `scripts/06_compute_cx_gene_matrix.py` |

| 07 | compute | BAM -> FPKM | `scripts/07_compute_bam_to_fpkm.py` |

| 08 | compute | Differential expression (NB/other) helper | `scripts/08_compute_nb_de.py` |

| 09 | compute | Merge expression + pairwise DE pipeline | `scripts/09_compute_merge_expr_pairwise_de.py` |

| 10 | plot | Plot whole profiles in 2D/3D | `scripts/10_plot_whole_profile_2d3d.py` |

| 11 | plot | Track vs expression heatmap | `scripts/11_plot_track_vs_expr_heatmap.py` |

| 12 | plot | Histone clustering pipeline from matrices | `scripts/12_plot_histone_cluster_pipeline.py` |

| 13 | plot | Methylation profile 2D/3D | `scripts/13_plot_methylation_profile_2d3d.py` |

| 14 | plot | Methylation vs expression heatmap | `scripts/14_plot_meth_vs_expr_heatmap.py` |

| 15 | plot | Gene tracks (2D/3D) | `scripts/15_plot_gene_tracks_2d3d.py` |

| 16 | plot | Gene circle plot | `scripts/16_plot_gene_circle_plot.py` |


## Step-by-step manuals


Below is a merged manual. For steps without a dedicated markdown in `docs/`, please refer to the script `--help` and in-script docstrings.


---

# Step 01 — GFF3/GTF -> BED + gene length table


# OmicsCanvas Step 1: GFF3 → BED6 + CDS/exon length table (English guide)

Script: `1trans_gff_to_bed_genes_length_EN.py`

---

## 1. Purpose
This script parses a **GFF3 annotation** and outputs:
1) **BED6** coordinates for `gene` or `mRNA/transcript` features.
2) A **length table** summarizing total **CDS** or **exon** length.
3) If you set `--bed-feature gene` but CDS/exon `Parent` points to transcripts:
   - The script automatically builds a transcript→gene mapping from `mRNA/transcript` records
     (`ID=transcript`, `Parent=gene`)
   - Then aggregates transcript lengths up to gene level.

---

## 2. Input requirements
- Input must be a standard **GFF3** (9 columns, TAB-delimited)
- Must contain:
  - `gene` or `mRNA/transcript` records (for BED)
  - `CDS` or `exon` records (for length)
- The attribute column should include:
  - `ID=` (or your custom `--id-key`)
  - `Parent=` (or your custom `--parent-key`)

**Multi-parent supported**: `Parent=a,b,c` will be split automatically.

---

## 3. Outputs
### 3.1 BED6 (default: `gene.bed`)
Format: `chrom  start  end  ID  score  strand`
- By default, coordinates are kept as-is from GFF3.
- With `--bed-zero-based`:
  - `start = start - 1`
  - `end` unchanged
  - This converts GFF3 (1-based, inclusive) into BED (0-based, half-open)

### 3.2 Length table (default: `feature_length.tsv`)
Two columns (no header): `ID\tlength`
- Default: CDS length
- Use `--length-feature exon` for exon length
- Use `--merge-overlaps` to merge overlaps before summing (recommended)

---

## 4. Parameters (same as `-h`)
### Required
- `-i/--gff3`: input GFF3 file

### Output
- `-o/--out-bed`: output BED6 path
- `-l/--out-length`: output length table path

### Options
- `--bed-feature {gene,mRNA,transcript}`: which feature to export as BED
- `--length-feature {CDS,exon}`: which feature to use for length
- `--id-key`: attribute key for feature ID (default: `ID`)
- `--parent-key`: attribute key for parent ID (default: `Parent`)
- `--bed-zero-based`: convert start to BED 0-based
- `--score`: BED score column value (default: 10)
- `--merge-overlaps`: merge overlaps before length summation
- `--min-match`: minimum BED-vs-length ID overlap, otherwise abort

> Tip: if you process only a subset of genes, reduce `--min-match` (e.g. 100 or 10).

---

## 5. Examples
### 5.1 BED IDs as transcripts (mRNA)
```bash
python 1trans_gff_to_bed_genes_length_EN.py \
  -i annotation.gff3 \
  --bed-feature mRNA \
  --length-feature CDS \
  --bed-zero-based \
  -o transcript.bed \
  -l transcript_cds_length.tsv \
  --merge-overlaps \
  --min-match 1000
```

### 5.2 BED IDs as genes, aggregate transcript lengths to genes
```bash
python 1trans_gff_to_bed_genes_length_EN.py \
  -i annotation.gff3 \
  --bed-feature gene \
  --length-feature CDS \
  --bed-zero-based \
  -o gene.bed \
  -l gene_cds_length.tsv \
  --merge-overlaps \
  --min-match 1000
```

### 5.3 Non-standard attribute keys
```bash
python 1trans_gff_to_bed_genes_length_EN.py \
  -i annotation.gff3 \
  --bed-feature gene \
  --id-key gene_id \
  --parent-key Parent \
  --length-feature exon \
  --bed-zero-based \
  -o gene.bed \
  -l gene_exon_length.tsv \
  --min-match 100
```

---

## 6. Troubleshooting
### "ID mismatch" / overlap too small
- Ensure `--bed-feature` matches the ID layer you want.
- Ensure CDS/exon `Parent` points to transcript IDs (common) or gene IDs.
- Double-check `--id-key` / `--parent-key`.
- Reduce `--min-match` if working with a subset.

### BED coordinate looks off
- Add `--bed-zero-based` if you need standard BED coordinates.



---

# Step 02 — Split CX report by context (CG/CHG/CHH)


# BamTrackSuite (CX Coverage Tools) User Guide (English)

This guide covers two utilities (Chinese/English scripts share identical logic; only comments/help language differs):

- `cx_context_split.en.py` / `cx_context_split.zh.py`: split Bismark `CX_report` into **CG/CHG/CHH** coverage files.
- `cx_replicate_merge.en.py` / `cx_replicate_merge.zh.py`: merge **2 or 3 replicates** by `(chrom, pos)` and sum columns **3/4**.

---

## 1. Requirements

### 1.1 Python
- Recommended: Python **3.8+**

### 1.2 Dependencies
- `cx_context_split.*.py`: requires `pandas`
- `cx_replicate_merge.*.py`: standard library only

Install:
```bash
pip install pandas
```

---

## 2. Tool 1: Split `CX_report` (CG/CHG/CHH)

### 2.1 What it does
Reads Bismark `CX_report` (TSV, no header) and writes three files:
- `<prefix>_CG.CX`
- `<prefix>_CHG.CX`
- `<prefix>_CHH.CX`

Each output line has 4 columns:
```
chrom    position    methylated_count    total_count
```
where `total_count = methylated + unmethylated`.

### 2.2 Input format
Typical Bismark `CX_report` columns:
1) chromosome
2) position (1-based)
3) strand
4) count_methylated
5) count_unmethylated
6) context (CG/CHG/CHH/...)

The script reads columns `[0,1,3,4,5]`, so the file must have at least 6 columns.

> Note: this script does **NOT** read `.gz` directly. Please decompress first:
```bash
zcat sample.CX_report.txt.gz > sample.CX_report.txt
```

### 2.3 Arguments
```bash
python cx_context_split.en.py -h
```
- `-i/--input-cx FILE` (required): input `CX_report` path.
- `-p/--prefix STR` (required): output prefix (sample name).
- `-d/--out-dir DIR` (default: `meth_data`): output directory.
- `--max-depth INT` (default: `300`): cap `total_count` to reduce extreme coverage impact; also enforces `methylated_count <= total_count`.
- `--chunksize INT` (default: `1000000`): pandas chunk size for large files.

### 2.4 Example
```bash
python cx_context_split.en.py \
  -i sample.CX_report.txt \
  -p sample \
  -d meth_data \
  --max-depth 300 \
  --chunksize 1000000
```

---

## 3. Tool 2: Merge replicates (2 or 3 inputs)

### 3.1 What it does
Merges 2 or 3 replicate files by `(chrom, pos)` and sums columns 3 and 4.

Input: at least 4 columns per line
```
chrom    pos    col3    col4
```
Output: exactly 4 columns
```
chrom    pos    sum_col3    sum_col4
```
Missing positions in some replicates are treated as 0.

### 3.2 Critical requirement: inputs must be sorted
Each input file must be sorted by `(chrom, pos)` ascending.

Recommended sorting (GNU sort):
```bash
sort -k1,1V -k2,2n input.CX > input.sorted.CX
sort -k1,1V -k2,2n -c input.sorted.CX
```

### 3.3 Arguments
```bash
python cx_replicate_merge.en.py -h
```
- `-o/--out FILE` (required): output path; `.gz` suffix enables gzip output.
- `inputs` (required): 2 or 3 input paths (plain text or `.gz`).

### 3.4 Examples
```bash
python cx_replicate_merge.en.py -o merged_CG.CX rep1_CG.CX rep2_CG.CX
python cx_replicate_merge.en.py -o merged_CG.CX.gz rep1_CG.CX.gz rep2_CG.CX.gz rep3_CG.CX.gz
```

---

## 4. Recommended workflow

Split each replicate first, then (optionally) sort and merge:
```bash
python cx_context_split.en.py -i rep1.CX_report.txt -p rep1 -d meth_data
python cx_context_split.en.py -i rep2.CX_report.txt -p rep2 -d meth_data

sort -k1,1V -k2,2n meth_data/rep1_CG.CX > meth_data/rep1_CG.sorted.CX
sort -k1,1V -k2,2n meth_data/rep2_CG.CX > meth_data/rep2_CG.sorted.CX

python cx_replicate_merge.en.py -o meth_data/merged_CG.CX \
  meth_data/rep1_CG.sorted.CX meth_data/rep2_CG.sorted.CX
```



---

# Step 03 — Merge CX replicates


**Phase:** prepare

Run:
```bash
python scripts/03_prepare_cx_replicate_merge.py --help
```



---

# Step 04 — Extract per-gene methylation table from CX


# OmicsCanvas – Extract Gene Methylation (CX) Records

This document describes how to use **`omicscanvas_extract_gene_methylation_EN.py`** to extract
per-base methylation (CX) records for **one gene** or a **gene list**, producing a TSV that can be
used by downstream OmicsCanvas plotting scripts.

---

## 1. What this script does

- Reads **gene/transcript coordinates** from a **GFF3** file.
- Reads **CX methylation files** (precomputed) in the format: `ch  po  me  al`.
- For each target gene, extracts records in `[start-distance, end+distance]`.
- Converts genomic positions to **strand-aware relative coordinates** (so all genes align 5'→3').

It **does not** compute methylation; it only subsets your existing CX files.

---

## 2. Input requirements

### 2.1 GFF3
The script reads these columns:
- `seqid` (chromosome)
- `type` (feature type, e.g. mRNA or gene)
- `start`, `end`
- `strand` (+ or -)
- attributes (9th column)

You must specify:
- `--feature-type` (default: `mRNA`)
- `--attr-key` to extract the target ID from attributes (default: `ID`)

### 2.2 CX files
Expected naming under `--cx-dir`:

- `<SRR>_<context><cx-suffix>`

Default `cx-suffix` is `.CX`, so for example:
- `SRR9321764_CHH.CX`

Expected file format (no header by default):

```
ch\tpo\tme\tal
```

- `po` should be integer genomic coordinate (commonly 1-based)
- `me`, `al` can be float/int (counts)

---

## 3. Relative coordinate convention (important)

Let the extracted window be:
- `po1 = start - distance`
- `po2 = end + distance`

Relative coordinate is:
- `+` strand: `rel_po = po - po1`
- `-` strand: `rel_po = po2 - po`

So **rel_po always increases along transcription direction (5'→3')**.

---

## 4. Output

For each combination of (SRR/prefix, gene, context), the script writes:

- `<outdir>/<prefix>__<gene>__<context>.tsv`

Columns:

- `name`: gene ID
- `po`: relative coordinate
- `me`: methylated counts (numerator)
- `al`: total counts (denominator)

---

## 5. Common usage

### 5.1 Extract one gene from one SRR (all contexts)

```bash
python omicscanvas_extract_gene_methylation_EN.py \
  --gff3 genome/annotation.gff3 \
  --cx-dir meth/meth_data \
  --srr SRR9321764 \
  --gene Potri.001G055900.5.v3.0 \
  --contexts CG,CHG,CHH \
  --distance 2000 \
  --outdir single_gene
```

### 5.2 Extract a gene list from multiple SRRs

```bash
python omicscanvas_extract_gene_methylation_EN.py \
  --gff3 genome/annotation.gff3 \
  --cx-dir meth/meth_data \
  --srr SRR1,SRR2,SRR3 \
  --prefix WT,mut1,mut2 \
  --gene-list genes.txt \
  --contexts CHH \
  --distance 1000 \
  --outdir out_single_gene \
  --chunksize 2000000
```

### 5.3 If your CX files are gzipped

```bash
python omicscanvas_extract_gene_methylation_EN.py \
  --gff3 genome/annotation.gff3 \
  --cx-dir meth/meth_data \
  --cx-suffix .CX.gz \
  --srr SRR9321764 \
  --gene YOUR_GENE \
  --contexts CHH
```

---

## 6. Performance notes

- The script groups target genes by chromosome.
- For each chromosome, it reads only the **union range** covering all target genes on that chromosome.
- Use `--chunksize` to reduce memory; set `--chunksize 0` if you have enough RAM and want faster IO.

---

## 7. Troubleshooting

1) **No target genes found**
- Verify `--feature-type` (gene vs mRNA)
- Verify `--attr-key` matches your attribute field

2) **Missing CX file warnings**
- Ensure files are named exactly `<SRR>_<context><cx-suffix>`

3) Output is empty
- Your CX file may not contain records in that interval
- Check whether `po` is 1-based vs 0-based in your CX file



---

# Step 05 — BAM -> gene matrix (TSS/gene/TES)


# OmicsCanvas Step 2: BAM → Gene‑centric Coverage Matrices

This step converts a **sorted + indexed** BAM file into three **gene‑centric coverage matrices** used by OmicsCanvas plotting and clustering scripts.

## What this script does

Given:

- A **sorted, indexed** BAM (`.bam` + `.bai`)
- A **BED6** gene annotation file (`chrom, start, end, feature_id, score, strand`)

It outputs three matrices:

1) **TSS window matrix** (TSS ± distance)
2) **TES window matrix** (TES ± distance)
3) **Scaled gene profile matrix** (upstream distance + scaled gene body + downstream distance)

Each matrix is **one feature per row** (usually a gene), with columns:

- `ID` (feature_id)
- `bin_1 ... bin_N` (binned mean coverage)

## Output naming scheme

By default (`--naming standard`), the script writes:

- `<prefix>_gene_profile_matrix.tsv`
- `<prefix>_tss_matrix.tsv`
- `<prefix>_tes_matrix.tsv`

Why this is more规范 than the old names:

- `start/end` are ambiguous; **TSS/TES** are standard biological terms.
- The "gene" matrix includes **upstream + scaled gene body + downstream**, so calling it only `gene_body` is slightly misleading.

For full backward compatibility with existing pipelines, use:

- `--naming legacy`

which writes:

- `<prefix>_gene_body_matrix.txt`
- `<prefix>_start_matrix.txt`
- `<prefix>_end_matrix.txt`

You can also override each suffix individually with `--suffix-tss / --suffix-gene / --suffix-tes`.

## Coordinate conventions (important)

- The script assumes **BED coordinates are 0‑based, half‑open** (standard BED): `[start, end)`.
- Coverage is computed on the BAM reference coordinates (0‑based indexing in Python slicing).

If your annotation is 1‑based (common in GFF/GTF), convert it to BED first (Step 1).

## Requirements

- Python ≥ 3.8
- `pysam`, `numpy`, `pandas`

The BAM must be:

- sorted: `samtools sort`
- indexed: `samtools index`

## Command‑line usage

### Basic (standard naming)

```bash
python omicscanvas_bam_to_gene_matrices_EN.py \
  --bam sample.sorted.bam \
  --bed genes.bed \
  --outdir caculate_matrix \
  --out-prefix B_H3K4me3 \
  --threads 12
```

This produces:

- `caculate_matrix/B_H3K4me3_gene_profile_matrix.tsv`
- `caculate_matrix/B_H3K4me3_tss_matrix.tsv`
- `caculate_matrix/B_H3K4me3_tes_matrix.tsv`

### Backward‑compatible output names (legacy)

```bash
python omicscanvas_bam_to_gene_matrices_EN.py \
  --bam sample.sorted.bam \
  --bed genes.bed \
  --out-prefix caculate_matrix/B_H3K4me3 \
  --naming legacy \
  --threads 12
```

### Custom suffixes

```bash
python omicscanvas_bam_to_gene_matrices_EN.py \
  --bam sample.sorted.bam \
  --bed genes.bed \
  --outdir caculate_matrix \
  --out-prefix B_H3K4me3 \
  --suffix-gene _gene_body_matrix.txt \
  --suffix-tss _start_matrix.txt \
  --suffix-tes _end_matrix.txt
```

## Parameter reference

### Required

- `--bam` : Sorted + indexed BAM.
- `--bed` : BED6 gene file: `chrom start end ID score strand`.
- `--out-prefix` : Output prefix (sample/mark name).

### Binning

- `--distance` (bp): The flank length around TSS/TES and upstream/downstream of the gene profile.
- `--bins-start`: #bins for **TSS ± distance** (total length = `2*distance`).
- `--bins-end`: #bins for **TES ± distance** (total length = `2*distance`).
- `--bins-gene-st`: #bins for **upstream distance** in the gene profile.
- `--bins-gene-body`: #bins to scale the **gene body** into (gene length can vary).
- `--bins-gene-en`: #bins for **downstream distance** in the gene profile.

**Divisibility constraints** (the script validates these):

- `(2*distance) % bins-start == 0`
- `(2*distance) % bins-end == 0`
- `distance % bins-gene-st == 0`
- `distance % bins-gene-en == 0`

### Normalization

Coverage is computed as:

```
raw_cov = A + C + G + T
norm_cov = raw_cov / mapped_reads / read_length / read_weight * norm_scale
```

- `--read-length`: read length used in the normalization.
- `--reads-type`: `single` or `paired` (paired uses weight=2).
- `--norm-scale`: multiplicative scale factor (default `1e8`).
- `--max-cov`: cap the normalized coverage to remove extreme outliers.

### Runtime and output

- `--threads`: number of worker processes.
- `--outdir`: output directory (optional).
- `--naming`: `standard` or `legacy` output naming.
- `--suffix-*`: override suffixes.

## Notes / pitfalls

1) Genes too close to chromosome ends may have incomplete windows and will be skipped (with a warning).
2) Ensure your BED chromosome names match BAM references (e.g., `Chr1` vs `1`).
3) For large genomes, this step can be heavy because it uses `pysam.count_coverage` per chromosome.



---

# Step 06 — Build CX_gene matrices/profiles from CX/genes


# CXGeneMatrix User Guide (English)

## 1. Overview
CXGeneMatrix converts a **single-context methylation CX file (CG/CHG/CHH)** into **feature-centric** methylation matrices and a global profile.

Typical use cases:
- Gene-centric metaplot matrices (upstream / gene body / downstream)
- TE-centric or any interval-centric profiles (as long as the BED format is valid)

Outputs:
1) per-feature bin tables (upstream / body / downstream)
2) a concatenated global profile (ratio = me/al per bin)

---

## 2. Requirements
- Python >= 3.8 (3.9+ recommended)
- Dependencies: numpy, pandas

Install:
```bash
pip install numpy pandas
```

---

## 3. Input formats

### 3.1 BED (>= 6 columns; strand is required)
Recommended 6 columns:
```
chrom   start   end     id      score   strand
Chr01   1000    2000    GeneA   0       +
Chr01   3000    3500    GeneB   0       -
```

Notes:
- The strand column must be `+` or `-`
- BED coordinates must be consistent with your CX file coordinates (BED is often 0-based; CX may be 1-based)
- You can pass a **TE BED** to generate TE-centric methylation matrices

---

### 3.2 CX (4 columns, no header)
```
ch      pos     me      al
Chr01   1050    3       10
Chr01   1051    0       8
Chr01   1052    1       6
```

Recommendation:
- Sort by (ch, pos). The script sorts after loading, but sorting a very large file in pandas can be expensive.

---

## 4. Arguments (args)

Run `python cx_gene_matrix.en.py -h` for a full help message.

### Required
- `-s / --sample`  
  Output prefix (sample/tissue name), e.g. `SRR9321764`, `S2`.
- `-b / --bed`  
  BED path (>= 6 columns, includes strand).

### Common
- `-c / --context`  
  Context: `CG` / `CHG` / `CHH` (default: CHH).
- `--cx-file`  
  Explicit CX path (ch, pos, me, al).
- `--cx-dir`  
  CX directory (default: `meth_data`). If `--cx-file` is omitted, CX is inferred as:  
  `<cx-dir>/<sample>_<context>.CX`
- `--distance`  
  Flanking extension in bp (default: 2000).
- `--bins-up` / `--bins-body` / `--bins-down`  
  Bin counts for upstream / body / downstream (default: 50 / 100 / 50).
- `--chrom-prefix`  
  Keep only chromosomes starting with this prefix (default: `Chr`). Disable by setting it to empty: `--chrom-prefix ""`.
- `-o / --out-dir`  
  Output directory (default: `CX_gene`).
- `--overwrite`  
  Replace existing outputs; otherwise the script fails if outputs already exist.

---

## 5. Outputs
Example: `sample=SRR9321764`, `context=CHH`
- `CX_gene/SRR9321764_CHH_upstream_bins50.tsv`
- `CX_gene/SRR9321764_CHH_body_bins100.tsv`
- `CX_gene/SRR9321764_CHH_downstream_bins50.tsv`
- `CX_gene/SRR9321764_CHH_profile.tsv`

Columns in the three matrix tables:
- `id`: feature ID from BED column 4
- `bin`: bin index starting from 1
- `me`: summed methylated count in the bin
- `al`: summed total coverage count in the bin

Profile:
- `ratio`: concatenated (upstream→body→downstream) `me/al` values

---

## 6. Examples

### 6.1 Gene-centric matrix
```bash
python cx_gene_matrix.en.py \
  -s SRR9321764 -c CHH \
  -b genes.bed \
  --cx-dir meth_data \
  --distance 2000 --bins-up 50 --bins-body 100 --bins-down 50 \
  -o CX_gene \
  --overwrite
```

### 6.2 TE-centric matrix (TE BED)
```bash
python cx_gene_matrix.en.py \
  -s SRR9321764 -c CHH \
  -b TE.bed \
  --cx-dir meth_data \
  --distance 2000 --bins-up 50 --bins-body 100 --bins-down 50 \
  -o CX_TE \
  --overwrite
```

---

## 7. Troubleshooting
1) **Missing strand in BED**
- Add the 6th column as `+/-`. If strand is unknown and orientation is not critical, use `+`.

2) **Empty outputs or many NaNs**
- Check chromosome naming and coordinate conventions; disable `--chrom-prefix` if needed.

3) **Memory issues for very large CX**
- Pre-sort CX with system tools, then run the script.



---

# Step 07 — BAM -> FPKM


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



---

# Step 08 — Differential expression (NB/other) helper


**Phase:** compute

Run:
```bash
python scripts/08_compute_nb_de.py --help
```



---

# Step 09 — Merge expression + pairwise DE pipeline


**Phase:** compute

Run:
```bash
python scripts/09_compute_merge_expr_pairwise_de.py --help
```



---

# Step 10 — Plot whole profiles in 2D/3D


# OmicsCanvas Step 3: Plot Whole Profiles (2D / 3D)

This step visualizes **meta‑profiles** (1D tracks) from the matrix files generated by the OmicsCanvas matrix generator (recommended) or the legacy generator.

- Script: `omicscanvas_plot_whole_profile_2d3d.py` (English)
- Chinese script: `omicscanvas_plot_whole_profile_2d3d_CN.py`

---

## 1. What this script does

Given one or more **matrix files** per sample (e.g., ChIP‑seq coverage, ATAC‑seq coverage, RNA coverage, methylation, etc.), the script:

1. Loads matrices by **sample prefix** and **region type** (TSS / gene / TES).
2. Applies optional gene ID filtering (default keeps only `.1` transcript IDs).
3. Computes a 1D profile per sample (mean/median across genes).
4. Plots profiles as either:
   - **2D**: one axis per panel; panels arranged in a grid (columns × panels).
   - **3D**: “stacked” axes (offset) inside each column with vertical connectors.

---

## 2. Input matrix naming (standard vs legacy)

This script supports two filename conventions controlled by `--naming`.

### Standard naming (recommended)
Default when `--naming standard`:

- **TSS** window: `*_tss_matrix.tsv`
- **gene profile**: `*_gene_profile_matrix.tsv`
- **TES** window: `*_tes_matrix.tsv`

### Legacy naming
When `--naming legacy`:

- **TSS** window: `*_start_matrix.txt`
- **gene profile**: `*_gene_body_matrix.txt`
- **TES** window: `*_end_matrix.txt`

You normally **do not** need to override suffixes manually. If you must, advanced options exist:

- `--suffix-tss` / `--suffix-start`
- `--suffix-gene`
- `--suffix-tes` / `--suffix-end`

---

## 3. The key arguments

### Required

- `--mode {2d,3d}`
  - `2d`: one axes per panel.
  - `3d`: stacked axes + vertical connectors (better for “floating layers”).

- `--matrix-dir DIR`
  Folder containing the matrix files.

- `--group STRING`
  **Panel layout + sample prefixes**. Separators:

  - `,` = multiple samples in the same panel (multiple lines)
  - `;` = multiple panels stacked in the same column
  - `|` = new column

  Each token is a **sample prefix** used to load: `<matrix-dir>/<prefix><suffix>`.

- `--names STRING`
  Legend names for each sample line. Must match the exact structure of `--group`.

- `--ylabels STRING`
  One y‑axis label per panel. Uses `;` and optional `|` (no comma split).

- `--out FILE`
  Output figure path (`.pdf`, `.png`, `.svg`, ...).

### Common scientific parameters

- `--gene-type {TSS,gene,TES}`
  Select which region to plot.

- `--distance INT`
  Flank size (bp) used only for axis labels (e.g., “±2000 bp”).

- `--bins-start INT`, `--bins-end INT`
  Total bins for TSS and TES windows.

- `--bins-gene-st INT`, `--bins-gene-body INT`, `--bins-gene-en INT`
  Bins for the three parts in the gene profile.

- `--stat {mean,median}`
  Statistic used to collapse genes → one profile.

- `--scale FLOAT`
  Multiply matrix values by this factor before profiling (useful when matrices store small decimals).

- `--index-filter REGEX`
  Filter which gene/transcript IDs to include (default keeps `.1`). Use `''` to disable.

### Appearance

- `--fig-x FLOAT`, `--fig-y FLOAT`
- `--line-lw FLOAT`
- `--cmap STR` (matplotlib colormap)
- `--line-colors "c1,c2,..."` (optional explicit colors)
- `--ylim "ymin,ymax;ymin,ymax;..."` per panel (optional)
- `--legend / --no-legend`, `--legend-loc {inside,outside}`, `--legend-anchor "x,y"`
- `--share-y {none,column,all}`

### Layout controls

These are useful when you want publication‑grade geometry.

- `--base-left`, `--base-top`
- `--col-width`, `--panel-height`
- `--col-gap`, `--panel-gap`
- `--x-off`, `--y-off` (3D stacking offsets)
- `--strict-layout` (raise if the grid would overflow)

---

## 4. JSON config to reduce repetitive arguments

You can put common arguments into a JSON file and load it with `--config`.

- Keys may use `-` or `_`.
- CLI options always override config values.

Example:

```bash
python omicscanvas_plot_whole_profile_2d3d.py \
  --config plot_profile_config.json \
  --mode 3d \
  --gene-type gene \
  --group "G_H3K4me3,B_H3K4me3,R_H3K4me3;G_H3K27me3,B_H3K27me3,R_H3K27me3" \
  --names "G,B,R;G,B,R" \
  --ylabels "H3K4me3;H3K27me3" \
  --out whole_profile_gene_3d.pdf
```

---

## 5. Worked examples

### Example A: 2D, one column, multiple panels

```bash
python omicscanvas_plot_whole_profile_2d3d.py \
  --mode 2d \
  --matrix-dir caculate_matrix \
  --gene-type gene \
  --group "G_Input,B_Input,R_Input;G_H3K4me3,B_H3K4me3,R_H3K4me3" \
  --names "G,B,R;G,B,R" \
  --ylabels "Input;H3K4me3" \
  --out whole_profile_gene_2d.png
```

### Example B: 3D, multiple columns

```bash
python omicscanvas_plot_whole_profile_2d3d.py \
  --mode 3d \
  --matrix-dir caculate_matrix \
  --gene-type TES \
  --group "G_H3K4me3,B_H3K4me3,R_H3K4me3|G_H3K27me3,B_H3K27me3,R_H3K27me3" \
  --names "G,B,R|G,B,R" \
  --ylabels "H3K4me3|H3K27me3" \
  --out whole_profile_tes_3d.pdf
```

### Example C: legacy matrices

```bash
python omicscanvas_plot_whole_profile_2d3d.py \
  --mode 2d \
  --naming legacy \
  --matrix-dir caculate_matrix \
  --gene-type gene \
  --group "LSH10_H3K4me3,WT_H3K4me3" \
  --names "LSH10,WT" \
  --ylabels "H3K4me3" \
  --out legacy_profile_gene.png
```

---

## 6. Troubleshooting

- **File not found**: check `--matrix-dir`, `--naming`, and your prefixes in `--group`.
- **Structure mismatch**: `--names` must mirror `--group` exactly.
- **Empty plots**: common reasons are `--index-filter` removing all genes, or matrices containing only NaNs.
- **Layout overflow**: adjust `--col-width`, `--panel-height`, `--col-gap`, `--panel-gap`, or disable `--strict-layout`.



---

# Step 11 — Track vs expression heatmap


# OmicsCanvas: Track vs Expression Heatmap (EN)

Script:
- `omicscanvas_histone_vs_expr_heatmap.py` (English help)
- `omicscanvas_histone_vs_expr_heatmap_CN.py` (Chinese help)

This step plots **ChIP-seq / ATAC-seq (or any binned track) matrices** as heatmaps after **binning genes by RNA expression**.
It is used to answer questions like:
- *How do different histone marks behave from TSS → gene body → TES as expression increases?*
- *Do highly expressed genes show stronger Pol II signal across the gene body?*

The script reads the **matrix files generated by Step2** (`omicscanvas_bam_to_gene_matrices.py`).

---

## 1) Required inputs

### 1.1 Matrix directory (Step2 output)
You pass a directory containing matrix files for each track/region.

**Recommended standardized suffixes (current OmicsCanvas convention):**
- `*_tss_matrix.tsv`
- `*_gene_profile_matrix.tsv`
- `*_tes_matrix.tsv`

Examples (prefix = `B_`):
- `B_H3K4me3_tss_matrix.tsv`
- `B_H3K4me3_gene_profile_matrix.tsv`
- `B_H3K4me3_tes_matrix.tsv`

If you still have old filenames, override with:
- `--suffix-tss`, `--suffix-gene`, `--suffix-tes`

### 1.2 Expression table (FPKM/TPM/count-derived)
Provide a TSV (recommended) where **row index is gene ID**.

Example (TSV):
```
GeneID	rep1	rep2	rep3
AT1G01010	5.2	4.9	5.1
AT1G01020	0	0	0
...
```

You select which columns are used to compute the mean expression with `--expr-cols`.

---

## 2) Track definition syntax (`--tracks`)

`--tracks` is **comma-separated**. Each item supports:

1) `NAME`
- Display name = `NAME`
- File ID = `<file-prefix>NAME`

2) `NAME:FILEID`
- Display name = `NAME`
- File ID = `<file-prefix>FILEID`

3) Replicates: `NAME:REP1+REP2+REP3`
- Each `REP*` becomes `<file-prefix>REP*`
- Replicates are averaged element-wise before plotting

Examples:

```bash
--tracks "ATAC,H3K4me1,H3K4me3" --file-prefix "B_"
```

```bash
--tracks "ATAC:ATAC_rep1+ATAC_rep2,H3K4me3" --file-prefix "B_"
```

---

## 3) Gene regions (`--gene-types`) and axis templates

`--gene-types` controls which region(s) are plotted:
- `TSS`  : around transcription start site
- `gene` : upstream + scaled gene body + downstream (full gene profile)
- `TES`  : around transcription end site

Example:
```bash
--gene-types TSS,gene,TES
```

The x-axis ticks/labels and vertical guide lines are derived from:
- `--distance`
- `--bins-start`, `--bins-end`
- `--bins-gene-st`, `--bins-gene-body`, `--bins-gene-en`

**Important:** these `--bins-*` must match the actual matrix column count.

---

## 4) Gene ID matching (expression ↔ matrices)

Different pipelines may use transcript IDs (e.g., `AT1G01010.1`) or gene IDs (`AT1G01010`).

Use:
- `--gene-id-mode exact`
  - IDs must match exactly
- `--gene-id-mode strip_dot`
  - `AT1G01010.1` → `AT1G01010`

If `strip_dot` leads to multiple isoforms per gene, `--collapse-isoform` controls how to merge:
- `mean` (default)
- `max`
- `first`

Optional filter:
- `--isoform-suffix .1`
  - only keep rows ending with `.1`

---

## 5) Expression binning

Genes are split into two parts:
- expression `<= --zero-threshold` → `--none-bins` groups
- expression `>  --zero-threshold` → `--exp-bins` groups

Each group produces **one heatmap row**.

Within each group, you choose the statistic across genes:
- `--stat mean` (default)
- `--stat median`

Sorting direction:
- default: low → high
- `--descending`: high → low

---

## 6) Color scale (vmin/vmax) and colormap

You can set color scaling in two ways:

### 6.1 Automatic scaling (recommended)
`--scale-mode`:
- `quantile` (default): use `--quantiles` (e.g., `0.01,0.99`)
- `diverging`: symmetric around 0, suitable for Z-score matrices
- `ratio`: 0..max, suitable for non-negative signal

Optional shortcut:
- `--range X`
  - diverging → `[-X, +X]`
  - others    → `[0, X]`

### 6.2 Manual scaling
- `--vmin` and `--vmax` override all automatic inference

Colormap:
- `--cmap RdBu_r` (diverging)
- `--cmap viridis` (sequential)

If `--cmap` is not set, the script chooses a reasonable default based on `--scale-mode`.

---

## 7) Outputs

For each (track × gene-type) combination, the script writes:

`<outdir>/<out-prefix>_<track>_<gene-type>_histone_vs_expr_heatmap.<ext>`

Where:
- `<ext>` is `pdf` or `png` controlled by `--out-format`

Plot controls:
- `--dpi` (PNG and PDF rasterization)
- `--ytick-step` (row tick spacing)
- `--no-border` (hide outer border)
- `--title` (optional title)

---

## 8) Example commands

### Example A: Standard naming, plot all regions
```bash
python omicscanvas_histone_vs_expr_heatmap.py \
  --matrix-dir caculate_matrix \
  --tracks "ATAC,H3K4me1,H3K4me3" \
  --file-prefix "B_" \
  --gene-types TSS,gene,TES \
  --expr FPKM_all.txt \
  --expr-cols 0,1,2 \
  --none-bins 10 \
  --exp-bins 90 \
  --distance 2000 \
  --bins-start 100 --bins-end 100 \
  --bins-gene-st 50 --bins-gene-body 100 --bins-gene-en 50 \
  --cmap RdBu_r \
  --scale-mode quantile --quantiles 0.01,0.99 \
  --out-format png \
  --outdir out_heatmap \
  --out-prefix B
```

### Example B: Old filenames (override suffix)
```bash
python omicscanvas_histone_vs_expr_heatmap.py \
  --matrix-dir caculate_matrix \
  --tracks "H3K4me3" --file-prefix "B_" \
  --gene-types gene \
  --suffix-gene _gene_body_matrix.txt \
  --expr FPKM_all.txt \
  --out-prefix legacy
```

---

## 9) Common pitfalls

1) **No overlap genes**
- Your matrix gene IDs and expression gene IDs differ (use `--gene-id-mode strip_dot` or adjust IDs).

2) **Bins mismatch**
- Your `--bins-*` settings do not match the matrix column count.

3) **Heatmap looks flat**
- Set a better scale: try `--scale-mode quantile --quantiles 0.05,0.95` or manually set `--vmin/--vmax`.



---

# Step 12 — Histone clustering pipeline from matrices


# OmicsCanvas – Cluster histone/track signals and link to expression (from Step2 matrices)

This script clusters genes using **ChIP-seq/ATAC-seq (or any track) gene-matrix outputs** from **OmicsCanvas Step2** (`omicscanvas_bam_to_gene_matrices.py`), then generates:

1. **Clustered heatmaps** (panel heatmaps + per-track heatmaps)
2. **Mean profile lines** for each cluster (per mark / per treatment)
3. **Expression distribution plots** (FPKM/TPM/CPM etc. in a TSV) per cluster

The workflow is designed for large matrices (many marks / treatments), so **panel heatmaps are saved as PNG** (PDF can be huge).

---

## Input files

### 1) Step2 matrix files
Generated by `omicscanvas_bam_to_gene_matrices.py`.

Each matrix is a TSV:
- **Index**: gene (or transcript) IDs
- **Columns**: bins, typically `bin_1 ... bin_N` (N depends on your bin settings)

**Standard naming (recommended / default):**
- TSS region: `<prefix>_tss_matrix.tsv`
- Gene profile region: `<prefix>_gene_profile_matrix.tsv`
- TES region: `<prefix>_tes_matrix.tsv`

These suffixes match your updated naming standard.

**Legacy naming (optional):**
- `<prefix>_start_matrix.txt`
- `<prefix>_gene_body_matrix.txt`
- `<prefix>_end_matrix.txt`

You can switch with `--naming legacy` or override suffixes with `--suffix-*`.

### 2) Expression table
A TSV where:
- **Column 1 (index_col=0)**: gene IDs
- Remaining columns: numeric expression columns (replicates)

You will map treatments to replicate columns using `--expr-cols`.

---

## How the grouping strings work

You control marks and treatments with **three aligned strings**:

### `--in-group` (prefix matrix groups)
- Split marks by `;`
- Within each mark, split treatments by `,`

Example with 2 marks and 3 treatments (G/B/R):

```
G_H3K4me3,B_H3K4me3,R_H3K4me3;G_H3K27me3,B_H3K27me3,R_H3K27me3
```

### `--in-ylabels` (mark labels)
Comma-separated mark names, one per `;` block in `--in-group`:

```
H3K4me3,H3K27me3
```

### `--in-names` (treatment labels)
- Split marks by `;`
- Within each mark, split treatments by `,`

Example:

```
G,B,R;G,B,R
```

---

## What the script does

### A) Load matrices
For each prefix in `--in-group`, the script loads:

`<matrix-dir>/<prefix><suffix>`

for regions selected in `--plot-regions` and for the region used in clustering (`--cluster-region`).

### B) ID normalization (matrix IDs ↔ expression IDs)
Step2 matrices sometimes use transcript IDs (e.g., `GeneX.1`) while expression uses gene IDs (`GeneX`).

Use `--id-mode`:
- `auto` (default): tries `none` → `drop_isoform` → `tx_to_gene`
- `drop_isoform`: remove isoform suffix like `.1`
- `tx_to_gene`: a common rule: `AT1G01010.1` → `AT1G01010`
- `regex`: custom mapping using `--id-regex` / `--id-repl`

If normalization creates duplicate rows, collapse with `--collapse-duplicates`.

### C) Build clustering feature matrix
Using only the region defined by `--cluster-region`:

1. For each track (mark × treatment), compute **zscore per gene across bins** (axis=1)
2. Concatenate all tracks into one large matrix (genes × features)
3. KMeans cluster (`--k`)

### D) Plot outputs
1) Panel heatmaps (2 PNGs each region: annotated + clean)
2) Per-track heatmaps (optional) + cluster boundaries
3) Mean profiles (per cluster) (PDF)
4) Expression boxplot per cluster (PDF)

---

## Typical command

```bash
python omicscanvas_histone_cluster_pipeline.py \
  --matrix-dir caculate_matrix \
  --in-group "G_H3K4me3,B_H3K4me3,R_H3K4me3;G_H3K27me3,B_H3K27me3,R_H3K27me3" \
  --in-ylabels "H3K4me3,H3K27me3" \
  --in-names "G,B,R;G,B,R" \
  --fpkm FPKM_all.txt \
  --expr-cols "G=0,1,2;B=3,4,5;R=6,7,8" \
  --cluster-region gene \
  --plot-regions TSS,gene,TES \
  --k 8 \
  --cmap RdBu_r \
  --zlim -2,2 \
  --out-prefix G_B_R_histone \
  --outdir histone_cluster_out
```

Notes:
- `--cluster-region gene` means clustering is computed from the **gene profile** region.
- `--plot-regions` controls which regions get heatmaps/profiles.

---

## Parameters (detailed)

### Input / naming
- `--matrix-dir` : folder containing Step2 matrices.
- `--naming` : `standard` (default) uses `_tss_matrix.tsv`, `_gene_profile_matrix.tsv`, `_tes_matrix.tsv`; `legacy` uses the old `*_start_matrix.txt` etc.
- `--suffix-tss`, `--suffix-gene`, `--suffix-tes` : override region suffixes (advanced).

### Grouping
- `--in-group` : prefixes grouped by mark and treatment.
- `--in-ylabels` : mark labels (comma-separated; align to blocks).
- `--in-names` : treatment labels per mark block.
- `--include-marks` / `--exclude-marks` : filter mark indices (0-based in the `;` order).

### Regions and bins
- `--cluster-region` : `TSS` / `gene` / `TES`.
- `--plot-regions` : comma-separated list of regions to plot.
- `--distance`, `--bins-start`, `--bins-gene-st`, `--bins-gene-body`, `--bins-gene-en`, `--bins-end` : only used to set tick marks / boundaries consistently with Step2.

### Normalization and clustering
- `--scale-factor` : multiply matrices by a constant after loading (kept for compatibility; default 10).
- `--zscore-axis` : typically 1 (per gene across bins).
- `--k` : number of clusters.
- `--random-state` : for reproducibility.

### Heatmap rendering
- `--cmap` : colormap (e.g., `RdBu_r`).
- `--zlim` : zscore range for heatmaps (`vmin,vmax`), e.g., `-2,2`.
- `--panel-cols`, `--panel-fig-x`, `--panel-fig-y`, `--panel-dpi` : panel heatmap layout.
- `--save-single` / `--no-single` : enable/disable per-track heatmaps.
- `--single-fig-x`, `--single-fig-y` : per-track heatmap size.

### Expression linking
- `--fpkm` : expression TSV.
- `--expr-cols` : map treatment name to replicate columns (0-based), e.g. `G=0,1,2;B=3,4,5;R=6,7,8`.
- `--expr-ylim` : y-limit for expression plot, e.g. `0,120`.

### ID alignment
- `--id-mode` : `auto|none|drop_isoform|tx_to_gene|regex`.
- `--id-regex`, `--id-repl` : only for `regex`.
- `--collapse-duplicates` : `mean|sum|first`.

### Output
- `--out-prefix` : prefix for output filenames.
- `--outdir` : output directory.

---

## Outputs
Inside `--outdir`:

- `figures/` : all plots
  - `*_panel_heatmaps_annot.png` : annotated panel heatmap
  - `*_panel_heatmaps_clean.png` : clean panel heatmap
  - `single_heatmaps/<region>/...png` : per-track heatmaps (if enabled)
  - `*_cluster_profiles.pdf` : per-cluster mean profiles
  - `*_cluster_expr_boxplot.pdf` : expression distribution per cluster

- `*_cluster_labels.tsv` : gene → cluster label
- `*_feature_matrix.tsv` : concatenated zscore feature matrix used for clustering

---

## Troubleshooting

- **No genes after merge**: ID mismatch. Try `--id-mode drop_isoform` or `--id-mode regex`.
- **Heatmaps look blank**: check `--zlim`, or confirm matrices are not all zeros.
- **Panel PNG too big**: reduce `--panel-fig-x/--panel-fig-y` or increase `--panel-cols`.



---

# Step 13 — Methylation profile 2D/3D


# OmicsCanvas — Methylation whole-profile plotting (2D / 3D)

This tool plots **methylation “whole-profile” curves** (typically CG/CHG/CHH) along a gene meta-profile, in either:
- **2D**: all contexts stacked vertically in one figure, each context showing multiple samples as lines
- **3D**: layered (pseudo-3D) stacking of contexts to visually separate layers

The script reads per-sample, per-context profile files and produces PDF figures.

---

## 1) Input files

### 1.1 File naming rule
For each sample and context, the script looks for:

```
<meth_dir>/<sample>_<context><whole_suffix>
```

Examples:
- `CX_gene/SRR9321764_CHH_profile.tsv`  (use `--meth-dir CX_gene --whole-suffix _profile.tsv`)
- `gene/LSH10_CG_whole_line.txt`        (use `--whole-suffix _whole_line.txt`)

### 1.2 File format
Each profile file should contain **one numeric column** (no header required). The script loads it as a single column named `methylation`.

The expected profile length is:

```
all_length = bins_st + bins_body + bins_en
```

If your profile length is different, set `--bins-st/--bins-body/--bins-en` to match how you generated the profile.

---

## 2) Outputs

Depending on `--mode`, the script writes:
- `<out_prefix>.2D.pdf`
- `<out_prefix>.3D.pdf`

Files are saved to the **current working directory**.

---

## 3) Common usage examples

### 3.1 Plot CHH profile for a single sample (2D)

```bash
python omicscanvas_methylation_profile_2d3d_EN.py \
  --meth-dir CX_gene \
  --contexts CHH \
  --samples SRR9321764 \
  --labels SRR9321764 \
  --whole-suffix _profile.tsv \
  --mode 2d \
  --out-prefix SRR9321764_CHH_profile
```

### 3.2 Plot CG/CHG/CHH for one sample (3D)

```bash
python omicscanvas_methylation_profile_2d3d_EN.py \
  --meth-dir CX_gene \
  --contexts CG,CHG,CHH \
  --samples SRR9321764 \
  --labels SRR9321764 \
  --whole-suffix _profile.tsv \
  --mode 3d \
  --out-prefix SRR9321764_CX_profile_3D
```

### 3.3 Plot multiple samples per context

```bash
python omicscanvas_methylation_profile_2d3d_EN.py \
  --meth-dir CX_gene \
  --contexts CHH \
  --samples SRR9321764,SRR9321765 \
  --labels WT,mut \
  --whole-suffix _profile.tsv \
  --mode both \
  --out-prefix CHH_WT_vs_mut
```

### 3.4 Provide manual line colors

```bash
python omicscanvas_methylation_profile_2d3d_EN.py \
  --meth-dir CX_gene \
  --contexts CHH \
  --samples SRR9321764,SRR9321765 \
  --labels WT,mut \
  --whole-suffix _profile.tsv \
  --line-colors '#1f77b4,#d62728' \
  --out-prefix CHH_color_demo
```

---

## 4) Parameter reference (public/common parameters)

### Input / sample definition
- `--meth-dir`  
  Directory that stores profile files.

- `--contexts`  
  Comma-separated contexts (e.g. `CG,CHG,CHH`). Each context becomes one layer/row.

- `--samples`  
  Comma-separated sample IDs used in the filename prefix `<sample>_<context>...`.

- `--labels`  
  Comma-separated legend labels for samples. Must match `--samples` length.

- `--whole-suffix`  
  Suffix appended after `<sample>_<context>`.  
  Examples: `_profile.tsv`, `_whole_line.txt`.

### Gene meta-profile layout (used for ticks/vertical lines)
- `--distance`  
  Only used for x-axis labeling (e.g. `-2000bp`/`+2000bp`). Not used for scaling.

- `--bins-st`  
  Number of bins in upstream (before TSS).

- `--bins-body`  
  Number of bins in gene body (TSS → TES).

- `--bins-en`  
  Number of bins in downstream (after TES).

### Figure & style
- `--fig-x`, `--fig-y`  
  Matplotlib figure size (inches).

- `--cmap`  
  Colormap name used to auto-assign line colors when `--line-colors` is not given.

- `--line-colors`  
  Comma-separated colors for samples, e.g. `'#1f77b4,#d62728'`. Colors cycle if fewer than samples.

- `--line-lw`  
  Line width.

- `--ylim`  
  Y-axis range(s). Either a single `min,max` for all contexts, or semicolon-separated per context.  
  Example: `0,0.3` or `0,0.3;0,0.2;0,0.1`.

- `--mode`  
  `2d`, `3d`, or `both`.

### 3D-only layout controls
- `--edge`  
  Outer margin (0–0.5 recommended). Larger = more whitespace.

- `--x-offset`, `--y-offset`  
  Layer offsets in percent (0–100). Controls pseudo-3D stacking.

### Output
- `--out-prefix`  
  Output file prefix.

- `--standard-distance-labels`  
  If set, TSS/TES labels are shown as `+distance ... -distance` style.

---

## 5) Troubleshooting

- **No output file appears**: outputs are written to the *current directory*; run `ls -lh *.pdf` after execution.
- **File not found**: confirm the naming rule `<sample>_<context><whole_suffix>` and set `--meth-dir` and `--whole-suffix`.
- **Wrong tick positions**: `--bins-st/--bins-body/--bins-en` must match the profile generation step.



---

# Step 14 — Methylation vs expression heatmap


# OmicsCanvas — Methylation vs Expression Heatmap

This tool integrates:
- **Methylation segment tables** (upstream / body / downstream)
- **Expression table (FPKM/TPM)**

and produces:
1) A **heatmap** showing methylation ratio along the merged profile after ranking genes by expression
2) A saved PDF figure for publication / review

---

## 1) Inputs

### 1.1 Methylation files (3 segments)
The script searches for three files:

```
<meth_dir>/<sample>_<cx><suffix-start>
<meth_dir>/<sample>_<cx><suffix-gene>
<meth_dir>/<sample>_<cx><suffix-end>
```

Example (your current naming):
- `CX_gene/SRR9321764_CHH_upstream_bins50.tsv`
- `CX_gene/SRR9321764_CHH_body_bins100.tsv`
- `CX_gene/SRR9321764_CHH_downstream_bins50.tsv`

Use:
- `--meth-dir CX_gene`
- `--suffix-start _upstream_bins50.tsv`
- `--suffix-gene  _body_bins100.tsv`
- `--suffix-end   _downstream_bins50.tsv`

#### Required columns
Each segment file must contain columns:

| column | meaning |
|---|---|
| `name` | gene ID |
| `po`   | bin index (can be any continuous integer range; the script will normalize by subtracting min) |
| `me`   | methylated counts |
| `al`   | total counts |

If the same `(name, po)` appears multiple times, the script sums `me` and `al` before computing ratios.

#### Segment length rules
- By default, the script requires the three segments to have the **same number of bins**.
- If your segments have different lengths (e.g., 50/100/50 is fine if each segment file itself has a consistent range), enable:

`--allow-unequal-segments`

The boundary positions (TSS/TES) will follow the inferred segment lengths.

### 1.2 Expression table
A TSV file where:
- **Row index (first column)** is gene ID
- Columns are replicates or samples

Use:
- `--fpkm FPKM.txt`
- `--fpkm-cols 0,1,2` (0-based column indices to average)

---

## 2) What the script does

1) Load upstream/body/downstream methylation tables
2) Normalize each segment `po` to start from 0 (subtract segment `po_min`)
3) Offset segments and concatenate them into one global coordinate
4) Compute methylation ratio per `(gene, global_po)`:

`ratio = me / al`

5) Load expression, compute mean across `--fpkm-cols`
6) Rank genes by expression:
   - genes with expression <= `--zero-threshold` go into `--none-bins` bins
   - genes with expression >  `--zero-threshold` go into `--exp-bins` bins
7) For each expression bin, aggregate methylation ratio along global bins
8) Plot a heatmap

---

## 3) Outputs

The figure is written to:

`<out_prefix>_<sample>_<cx>_meth_vs_expr_heatmap.pdf`

Saved in the current working directory.

---

## 4) Example command (your current naming)

```bash
python omicscanvas_meth_vs_expr_heatmap_EN.py \
  --sample SRR9321764 \
  --cx CHH \
  --meth-dir CX_gene \
  --suffix-start _upstream_bins50.tsv \
  --suffix-gene  _body_bins100.tsv \
  --suffix-end   _downstream_bins50.tsv \
  --fpkm FPKM.txt \
  --fpkm-cols 0,1,2 \
  --none-bins 10 \
  --exp-bins 90 \
  --distance 2000 \
  --scale-mode ratio \
  --cmap RdBu_r \
  --vmax 0.25 \
  --out-prefix SRR9321764
```

---

## 5) Parameter reference (common ones)

- `--sample` : sample ID for filename prefix.
- `--cx` : context (CG/CHG/CHH).
- `--meth-dir` : directory for methylation segment tables.
- `--suffix-start/--suffix-gene/--suffix-end` : file suffix for each segment.
- `--fpkm` : expression table.
- `--fpkm-cols` : 0-based indices to average.
- `--none-bins` : number of bins for genes with expression <= threshold.
- `--exp-bins` : number of bins for genes with expression > threshold.
- `--zero-threshold` : threshold (default 0.0).
- `--distance` : label value shown at profile edges (e.g. -2000bp / +2000bp).
- `--scale-mode` :
  - `ratio` (default): use raw ratio, typically [0, 1]
  - `diverging`: symmetric scale around 0 (useful after transforming)
  - `quantile`: auto-range from quantiles
- `--vmin/--vmax` : manually set color scale bounds.
- `--range` : convenience bound (if given, sets both vmin/vmax depending on mode).
- `--cmap` : colormap.
- `--out-prefix` : output prefix.

---

## 6) Troubleshooting

- **Heatmap is empty/blank**: most often caused by no overlap between methylation gene IDs and expression gene IDs.
- **Division by zero**: bins with `al=0` become NaN; you can fill with `--fillna 0`.
- **Segment bins error**: if segment lengths differ, enable `--allow-unequal-segments`.



---

# Step 15 — Gene tracks (2D/3D)


# OmicsCanvas — Single-gene Track Plot (2D/3D)

This script draws **stacked linear tracks** for a single gene: gene model (exons/introns) + multiple **BAM-derived coverage tracks** + optional **methylation tracks**. It supports **2D stacked** and **3D layered** rendering.

- Script (EN): `omicscanvas_gene_tracks_2d3d_EN.py`
- Script (CN): `omicscanvas_gene_tracks_2d3d_CN.py`

---

## 1) Inputs

### A) Genome annotation (GFF3)
Provide a GFF3 file for drawing the gene structure and for locating the target gene.

- `--gff3`: path to `*.gff3` (or `*.gff.gz` if your Python can read it)
- `--gene`: gene ID to plot (e.g., `Potri.001G055900.5.v3.0`)

If your GFF3 uses non-standard feature names, use:
- `--gene-keywords`: keywords/regex used to detect gene/transcript features (default: `gene|mRNA|transcript`)
- `--exon-keywords`: keywords/regex used to detect exon/CDS features (default: `exon|CDS`)

### B) BAM coverage tracks
Coverage is computed by reading BAM files around the gene region.

- `--bam-dir`: directory containing BAMs (default: `./bam`)
- `--bam-spec`: **track grouping spec** (see below)
- `--name-spec`: labels for each track in `--bam-spec` (same shape)

**Track grouping spec format**
- Use `;` to separate **groups** (rows / layers)
- Use `,` to separate **tracks** inside a group

Example:

```txt
--bam-spec  "H3K27me3.bam,H3K36me3.bam;ATAC.bam;RNA.bam" \
--name-spec "H3K27me3,H3K36me3;ATAC;RNA"
```

### C) Optional methylation tracks (prefix-based)
If provided, the script will load methylation TSVs by **prefix**.

- `--meth-dir`: folder containing methylation TSVs
- `--meth-spec`: same grouping syntax as BAM, but elements are **prefixes** (not full filenames)

The script expects per-prefix files in the form:
- `<prefix>__CG.tsv`
- `<prefix>__CHG.tsv`
- `<prefix>__CHH.tsv`

Example:

```txt
--meth-dir CX_gene \
--meth-spec "SRR9321764_CHH"  
```

> If your methylation TSV uses a different naming convention, change `load_meth_prefix()` inside the script.

---

## 2) Output

- `--out`: output figure path (recommended: `.png` for quick viewing or `.pdf` for vector)

Notes:
- **Large multi-track figures** can be huge as PDF. Prefer `png` for very large layouts.
- Use `--fig-x/--fig-y` to control the canvas size.

---

## 3) Key parameters

### Rendering
- `--mode {2d,3d}`: choose 2D stacked or 3D layered view
- `--fig-x`, `--fig-y`: figure size in inches
- `--cmap`: colormap for track heatmaps (default: `RdBu_r`)

### Genomic window / bins
- `--distance`: flank length around the gene (bp). For example `2000` means `-2000bp .. +2000bp` around TSS/TES.

### Coverage normalization
- `--reads-length`: approximate fragment length for extension (default: 200)
- `--seq-type {single,paired}`: how to interpret reads (default: `paired`)

### Y-axis controls (2D/3D)
- `--ylim "min,max"`: global y-limit for all tracks (e.g. `0,50`)
- `--group-ylims`: per-group y-limits using `;` (e.g. `0,10;0,50;0,5`)
- `--share-group-ylim`: share one y-limit per group

### Methylation options
- `--meth-layout {combined,separate}`: combined = CG/CHG/CHH merged; separate = 3 subtracks
- `--meth-ylabel`: y-axis label for methylation (default: `Methylation`)
- `--meth-hide-prefix`: hide prefix text in methylation labels
- `--meth-drop-me-zero`: drop sites with `me==0` before plotting (recommended)

### 3D layout (only used in `--mode 3d`)
These control the 3D “panel” placement.
- `--base-left`, `--top-start`
- `--track-width`, `--track-height`
- `--x-off`, `--y-off`
- `--layer-gap`

If you don’t want to tune these, start with defaults and only adjust figure size.

---

## 4) Examples

### Example 1 — 2D stacked track plot
```bash
python omicscanvas_gene_tracks_2d3d_EN.py \
  --mode 2d \
  --gff3 genome.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "H3K27me3.bam,H3K36me3.bam;ATAC.bam;RNA.bam" \
  --name-spec "H3K27me3,H3K36me3;ATAC;RNA" \
  --meth-dir CX_gene \
  --meth-spec "SRR9321764" \
  --meth-layout combined \
  --fig-x 12 --fig-y 6 \
  --out gene_tracks_2d.png
```

### Example 2 — 3D layered view
```bash
python omicscanvas_gene_tracks_2d3d_EN.py \
  --mode 3d \
  --gff3 genome.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "H3K27me3,H3K36me3;ATAC;RNA" \
  --name-spec "H3K27me3,H3K36me3;ATAC;RNA" \
  --fig-x 14 --fig-y 7 \
  --out gene_tracks_3d.png
```

---

## 5) Troubleshooting

- **No reads / empty tracks**: check BAM index (`.bai`) exists and region names match the GFF3 (chr naming).
- **Gene not found**: adjust `--gene-keywords/--exon-keywords`, or pass the transcript ID if your file is transcript-centric.
- **Layout overlaps**: increase `--fig-x/--fig-y`, then tune `--track-width/--track-height` (3D) or `--ylim`.



---

# Step 16 — Gene circle plot


# OmicsCanvas — Single-gene Circle Plot

This script renders a **circular** view of a single gene, with optional **BAM coverage rings** and **methylation rings**.

- Script (EN): `omicscanvas_gene_circle_plot_EN.py`
- Script (CN): `omicscanvas_gene_circle_plot_CN.py`

---

## 1) Inputs

### A) GFF3 + Gene ID
- `--gff3`: genome annotation (GFF3)
- `--gene`: gene/transcript ID
- `--distance`: flank length (bp)

### B) BAM rings
- `--bam-dir`: directory containing BAMs
- `--bam-spec`: `;` groups and `,` tracks (same syntax as the 2D/3D script)
- `--name-spec`: labels aligned to `--bam-spec`

### C) Methylation rings (optional)
- `--meth-dir`: directory of methylation TSVs
- `--meth-spec`: prefix list (same grouping syntax)

Expected prefix files:
- `<prefix>__CG.tsv`, `<prefix>__CHG.tsv`, `<prefix>__CHH.tsv`

---

## 2) Output
- `--out`: output path (recommend `.png` for convenience)

---

## 3) Key parameters

### Figure
- `--fig-x`, `--fig-y`: figure size in inches

### BAM processing
- `--reads-length`: fragment length extension (default 200)
- `--seq-type {single,paired}`

### Circle layout
- `--circle-bam-step`: down-sampling step for BAM points on the ring (smaller = denser)
- `--circle-ring-start`: inner radius to start drawing rings
- `--circle-ring-width`: width of each ring
- `--circle-ring-gap`: gap between rings

### Style
- `--circle-frame / --no-circle-frame`: draw outer frame
- `--circle-frame-lw`: frame line width
- `--circle-margin`: extra margin
- `--circle-lim`: manual color/height scaling for rings

### Methylation
- `--circle-meth-layout {combined,separate}`: combined = merge CG/CHG/CHH; separate = three rings
- `--meth-drop-me-zero`: drop `me==0` sites

---

## 4) Example
```bash
python omicscanvas_gene_circle_plot_EN.py \
  --gff3 genome.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "H3K27me3.bam,H3K36me3.bam;ATAC.bam" \
  --name-spec "H3K27me3,H3K36me3;ATAC" \
  --meth-dir CX_gene \
  --meth-spec "SRR9321764" \
  --circle-meth-layout combined \
  --fig-x 10 --fig-y 10 \
  --out gene_circle.png
```

---

## 5) Troubleshooting
- If the circle looks too dense: increase `--circle-bam-step` or increase figure size.
- If rings are clipped: increase `--circle-margin` or disable frame.
