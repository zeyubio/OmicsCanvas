# 02_prepare_cx_context_split.py

OmicsCanvas / BamTrackSuite utility — split Bismark **CX_report** into **CG / CHG / CHH** context files.

This script reads a Bismark `CX_report` (TSV) and writes **three coverage files**:

- `<prefix>_CG.CX`
- `<prefix>_CHG.CX`
- `<prefix>_CHH.CX`

Each output line has **4 columns**:

```
chrom    position    methylated_count    total_count
```

Where `total_count = methylated_count + unmethylated_count` (optionally depth-capped). fileciteturn55file7L47-L75

---

## 1) Requirements

- Python ≥ 3.8
- `pandas` fileciteturn55file7L48-L56

Install:
```bash
pip install pandas
```

---

## 2) Input format (CX_report)

Typical Bismark `CX_report` columns (no header):

1. chrom
2. pos (1-based)
3. strand
4. count_methylated
5. count_unmethylated
6. context (CG/CHG/CHH/...)

This script reads **columns [0,1,3,4,5]** (chrom, pos, me, un, ctx), so the file must have **at least 6 columns**. fileciteturn55file7L47-L56

### About `.gz`
The script relies on `pandas.read_csv()`. In most environments pandas can read `*.gz` transparently; if your environment cannot, decompress first (e.g., `zcat input.gz > input`). (Your older manual currently says “not supported”; that statement is outdated.) fileciteturn55file2L1-L4

---

## 3) Output behavior

- Output directory is created if missing (`--out-dir`). fileciteturn55file7L20-L26
- Existing output files are **deleted** before writing to avoid accidental appends. fileciteturn55file7L27-L31
- Rows with `total_count == 0` are dropped. fileciteturn55file7L62-L66
- If `--max-depth > 0`, the script caps `total_count` to `max_depth` and enforces `methylated_count <= total_count`. fileciteturn55file7L66-L70

---

## 4) Parameters (checked against the script)

Required:
- `-i / --input-cx FILE` : input `CX_report` path
- `-p / --prefix STR` : output prefix (sample name)

Optional:
- `-d / --out-dir DIR` (default: `meth_data`) : output directory
- `--max-depth INT` (default: `300`) : depth cap for `total_count`
  - Set to `0` to disable capping (because the script only caps when `max_depth > 0`). fileciteturn55file7L66-L70
- `--chunksize INT` (default: `1000000`) : pandas chunk size for streaming large files fileciteturn55file6L3-L9

---

## 5) Examples

### 5.1 Basic usage
```bash
python 02_prepare_cx_context_split.py \
  -i sample.CX_report.txt \
  -p sample \
  -d meth_data \
  --max-depth 300 \
  --chunksize 1000000
```

Outputs:
- `meth_data/sample_CG.CX`
- `meth_data/sample_CHG.CX`
- `meth_data/sample_CHH.CX` fileciteturn55file6L11-L15

### 5.2 Disable depth capping
```bash
python 02_prepare_cx_context_split.py \
  -i sample.CX_report.txt \
  -p sample \
  --max-depth 0
```

---

## 6) Troubleshooting

### “Input file not found”
Check the `-i/--input-cx` path. The script errors if the file does not exist. fileciteturn55file7L42-L44

### Output is empty
- Confirm the input really contains contexts `CG/CHG/CHH` in column 6.
- If your input uses lowercase contexts, normalize upstream (this script matches `CG/CHG/CHH` exactly). fileciteturn55file7L71-L74

---

## 7) Help
```bash
python 02_prepare_cx_context_split.py -h
```
