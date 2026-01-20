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
