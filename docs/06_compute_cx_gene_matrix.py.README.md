# 06_compute_cx_gene_matrix.py — CXGeneMatrix (feature-centric methylation matrices)

This script converts a **single-context methylation CX** file (CG/CHG/CHH) into:

1) **Per-feature methylation matrices** (upstream / body / downstream bin tables)  
2) A **global concatenated profile** (`ratio = me/al`) across bins (upstream→body→downstream)

It works with **genes** (BED6) and also any other interval set (e.g., TE BED), as long as the BED has strand. fileciteturn61file0L7-L23

---

## 1) Requirements

- Python ≥ 3.8
- `numpy`, `pandas` fileciteturn61file0L46-L47

```bash
pip install numpy pandas
```

---

## 2) Inputs

### 2.1 BED (required)

BED must have **at least 6 columns**:

```
chrom  start  end  id  score  strand
```

- `strand` must be `+` or `-` (required). fileciteturn61file5L24-L35

Coordinates:
- The script assumes your BED coordinates are consistent with CX coordinates. The code does **not** auto-convert between 0-based and 1-based. fileciteturn61file0L20-L23

### 2.2 CX (required)

CX must be 4 columns (no header):

```
ch  pos  me  al
```

- `pos` is integer; `me/al` are numeric. fileciteturn61file5L44-L55

Default CX path rule:
- If `--cx-file` is not provided, the script reads:  
  `<cx-dir>/<sample>_<context>.CX` fileciteturn61file2L26-L27

If your CX file is named differently (e.g., `.CX.gz`), pass the full path via `--cx-file`. fileciteturn61file7L42-L53

---

## 3) What the script outputs

For `sample=SRR9321764`, `context=CHH`, default bins (50/100/50):

- `<out-dir>/SRR9321764_CHH_upstream_bins50.tsv`
- `<out-dir>/SRR9321764_CHH_body_bins100.tsv`
- `<out-dir>/SRR9321764_CHH_downstream_bins50.tsv`
- `<out-dir>/SRR9321764_CHH_profile.tsv` fileciteturn61file0L31-L35 fileciteturn61file2L32-L37

Matrix tables contain:
- `id`: feature ID from BED column 4
- `bin`: bin index
- `me`: summed methylated counts in the bin
- `al`: summed total counts in the bin

Profile contains:
- `ratio`: a single column of length `bins_up + bins_body + bins_down`, produced by concatenating ratios from the three segments. fileciteturn61file2L1-L20

Overwrite behavior:
- If any output file exists and `--overwrite` is **not** set, the script fails with an error.
- With `--overwrite`, existing outputs are removed first. fileciteturn61file5L5-L13

---

## 4) Parameters (checked)

Required:
- `-s/--sample` : sample name used in output prefix fileciteturn61file7L24-L27
- `-b/--bed` : BED6 input fileciteturn61file7L33-L40

Common:
- `-c/--context {CG,CHG,CHH}` (default: CHH) fileciteturn61file7L28-L31
- `--cx-file FILE` : explicit CX path fileciteturn61file7L42-L48
- `--cx-dir DIR` (default: `meth_data`) : used to infer default CX path fileciteturn61file7L50-L53

Binning:
- `--distance` (default: 2000) fileciteturn61file7L55-L59
- `--bins-up` (default: 50) fileciteturn61file7L60-L63
- `--bins-body` (default: 100) fileciteturn61file7L64-L67
- `--bins-down` (default: 50) fileciteturn61file7L68-L71

Chromosome filter:
- `--chrom-prefix` (default: `Chr`)  
  Keep only records whose chromosome starts with this prefix; set empty to disable. fileciteturn61file7L73-L77

Output:
- `-o/--out-dir` (default: `CX_gene`) fileciteturn61file7L79-L83
- `--overwrite` : overwrite existing outputs fileciteturn61file7L84-L87

---

## 5) Examples (name unified)

### 5.1 Gene-centric matrices (default CX naming)
```bash
python 06_compute_cx_gene_matrix.py \
  -s SRR9321764 -c CHH \
  -b genes.bed \
  --cx-dir meth_data \
  --distance 2000 --bins-up 50 --bins-body 100 --bins-down 50 \
  -o CX_gene \
  --overwrite
```

### 5.2 TE-centric matrices (use TE BED)
```bash
python 06_compute_cx_gene_matrix.py \
  -s SRR9321764 -c CHH \
  -b TE.bed \
  --cx-dir meth_data \
  --distance 2000 --bins-up 50 --bins-body 100 --bins-down 50 \
  -o CX_TE \
  --overwrite
```

### 5.3 Non-standard CX filename (use --cx-file)
```bash
python 06_compute_cx_gene_matrix.py \
  -s SRR9321764 -c CHH \
  -b genes.bed \
  --cx-file meth_data/SRR9321764_CHH.CX.gz \
  -o CX_gene \
  --overwrite
```

---

## 6) Troubleshooting

- **“BED must have at least 6 columns”**: your BED is missing strand; add the 6th column `+/-`. fileciteturn61file5L24-L35
- **“No chromosomes left after filtering”**: disable `--chrom-prefix` or change it to match your chromosome names. fileciteturn61file2L50-L52
- **Many NaNs in profile**: `al==0` for those bins; verify coverage and coordinate consistency. fileciteturn61file2L8-L13

---

## 7) Help
```bash
python 06_compute_cx_gene_matrix.py -h
```

(Argparse `prog` is `CXGeneMatrix`.) fileciteturn61file7L18-L22
