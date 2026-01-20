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

