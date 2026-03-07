# OmicsCanvas Step 2: BAM → Gene-centric Coverage Matrices

Repo script: `05_compute_bam_to_gene_matrices.py`  
(Argparse `prog`: `omicscanvas_bam_to_gene_matrices`)

This step converts a **sorted + indexed** BAM file into three **gene-centric coverage matrices** used by OmicsCanvas plotting and clustering tools.

---

## 1) What this script does

**Inputs**
- A **sorted, indexed** BAM (`.bam` + `.bai`)
- A **BED6** annotation file: `chrom  start  end  ID  score  strand`

**Outputs** (one feature per row; `ID` + binned mean coverage)
- **TSS window** matrix (TSS ± distance)
- **Scaled gene profile** matrix (upstream distance + scaled gene body + downstream distance)
- **TES window** matrix (TES ± distance)

**Strand handling**
- For negative-strand features, bins are reversed so all outputs are aligned 5′→3′.

**Normalization**
Coverage is computed by `pysam.count_coverage` per chromosome:

```
raw_cov  = A + C + G + T
norm_cov = raw_cov / mapped_reads / read_length / lib_factor * norm_scale
lib_factor = 1 (single-end) or 2 (paired-end)
```

Values are capped at `--max-cov` (winsorization).

---

## 2) Coordinate conventions (important)

- BED coordinates are assumed **0-based, half-open**: `[start, end)`.
- If your annotation is 1-based (GFF/GTF style), convert to BED first (Step 1).

> Tip: chromosome names in BED **must** match BAM references (e.g., `Chr01` vs `1`), otherwise those chromosomes are skipped.

---

## 3) Output naming

### 3.1 Standard naming (default)
With `--naming standard` (default), outputs are:

- `<prefix>_tss_matrix.tsv`
- `<prefix>_gene_profile_matrix.tsv`
- `<prefix>_tes_matrix.tsv`

### 3.2 Legacy naming
With `--naming legacy`, outputs are:

- `<prefix>_start_matrix.txt`
- `<prefix>_gene_body_matrix.txt`
- `<prefix>_end_matrix.txt`

### 3.3 Custom suffixes
You can override each suffix regardless of naming mode:

- `--suffix-tss`
- `--suffix-gene`
- `--suffix-tes`

### 3.4 `--outdir` vs `--out-prefix`
- If `--out-prefix` contains a path separator (`/`), it is used **as-is** and `--outdir` is ignored.
- Otherwise, if `--outdir` is set, outputs are written as `outdir/<prefix>_*`.

---

## 4) Quick start

### 4.1 Basic run (standard naming)
```bash
python 05_compute_bam_to_gene_matrices.py \
  --bam sample.sorted.bam \
  --bed genes.bed \
  --outdir caculate_matrix \
  --out-prefix B_H3K4me3 \
  --threads 12
```

Produces:
- `caculate_matrix/B_H3K4me3_tss_matrix.tsv`
- `caculate_matrix/B_H3K4me3_gene_profile_matrix.tsv`
- `caculate_matrix/B_H3K4me3_tes_matrix.tsv`

### 4.2 Legacy output names (compat mode)
```bash
python 05_compute_bam_to_gene_matrices.py \
  --bam sample.sorted.bam \
  --bed genes.bed \
  --out-prefix caculate_matrix/B_H3K4me3 \
  --naming legacy \
  --threads 12
```

### 4.3 Custom suffixes
```bash
python 05_compute_bam_to_gene_matrices.py \
  --bam sample.sorted.bam \
  --bed genes.bed \
  --outdir caculate_matrix \
  --out-prefix B_H3K4me3 \
  --suffix-gene _gene_body_matrix.txt \
  --suffix-tss _start_matrix.txt \
  --suffix-tes _end_matrix.txt
```

---

## 5) Parameters (checked against the script)

### 5.1 Optional JSON config
- `--config CONFIG.json`  
  Optional JSON config file. CLI arguments override config values.

Defaults are taken from config if present; otherwise hard defaults are used:
- `distance=2000`
- `tss_bins=100`, `tes_bins=100`
- `gene_up_bins=50`, `gene_body_bins=100`, `gene_down_bins=50`
- `threads=5`

Example config (minimal):
```json
{
  "genomic_bins": {
    "distance": 2000,
    "tss_bins": 100,
    "tes_bins": 100,
    "gene_body": {
      "upstream_bins": 50,
      "body_bins": 100,
      "downstream_bins": 50
    }
  },
  "threads": 12,
  "paths": {
    "matrix_dir": "caculate_matrix"
  }
}
```

### 5.2 Required inputs
- `-b/--bam BAM` (required)  
  Sorted and indexed BAM (`.bam` + `.bai`).
- `-g/--bed BED` (required)  
  BED6: `chrom start end ID score strand(+/-)`.
- `-o/--out-prefix STR` (required)  
  Output prefix (sample/track name).

### 5.3 Binning / region definition
- `--distance INT`  
  Flank length (bp) used for TSS/TES windows and upstream/downstream of gene profile.

**Preferred bin parameter names (recommended)**
- `--tss-bins INT` : number of bins for **TSS ± distance** (total length `2*distance`)
- `--tes-bins INT` : number of bins for **TES ± distance** (total length `2*distance`)
- `--gene-up-bins INT` : number of bins for **upstream distance** in gene profile
- `--gene-body-bins INT` : number of bins for **scaled gene body**
- `--gene-down-bins INT` : number of bins for **downstream distance** in gene profile

**Deprecated (still accepted for backward compatibility; hidden in `-h`)**
- `--bins-start` → `--tss-bins`
- `--bins-end` → `--tes-bins`
- `--bins-gene-st` → `--gene-up-bins`
- `--bins-gene-body` → `--gene-body-bins`
- `--bins-gene-en` → `--gene-down-bins`
- `--gene-upstream-bins`, `--gene-downstream-bins` → same targets

**Divisibility constraints** (validated)
- `(2*distance) % tss_bins == 0`
- `(2*distance) % tes_bins == 0`
- `distance % gene_up_bins == 0`
- `distance % gene_down_bins == 0`

### 5.4 Normalization
- `--read-length INT` (default: `150`)  
  Read length used in normalization.
- `--reads-type {single,paired}` (default: `paired`)  
  Library type used in normalization (`paired` uses factor=2).
- `--norm-scale FLOAT` (default: `1e8`)  
  Multiplicative scale factor after normalization.
- `--max-cov FLOAT` (default: `500.0`)  
  Cap normalized coverage values at this threshold.

### 5.5 Runtime
- `-t/--threads INT`  
  Number of worker processes. (Default from config or `5`.)

### 5.6 Output control
- `--outdir DIR`  
  Output directory (optional; default may come from config `paths.matrix_dir`).
- `--naming {standard,legacy}` (default: `standard`)  
  Output filename scheme.
- `--suffix-tss`, `--suffix-gene`, `--suffix-tes`  
  Override suffixes.

---

## 6) Notes / pitfalls

1. **Genes near chromosome ends** may have incomplete windows and will be skipped silently (no row written).
2. This step can be heavy for large genomes because it runs `count_coverage` per chromosome.
3. Ensure BAM is truly indexed and has mapped reads; otherwise the script errors.

---

## 7) Help
```bash
python 05_compute_bam_to_gene_matrices.py -h
```
