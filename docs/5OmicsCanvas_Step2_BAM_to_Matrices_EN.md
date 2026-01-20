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
