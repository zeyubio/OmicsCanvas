# 03_prepare_cx_replicate_merge.py

OmicsCanvas / BamTrackSuite utility — **merge 2–3 replicate CX coverage files** by `(chrom, pos)`.

This script is designed to merge the **CG/CHG/CHH `.CX` files** produced by the context splitter (`02_prepare_cx_context_split.py`), but it works for any tab-delimited file whose first 4 columns are:

```
chrom   pos   col3   col4
```

It merges by key `(chrom, pos)` and outputs:

```
chrom   pos   sum_col3   sum_col4
```

- Lines starting with `#` and blank lines are ignored.
- Files may be plain text or gzip-compressed (`.gz`).
- If the output path ends with `.gz`, the output is gzip-compressed. fileciteturn56file0L1-L18

---

## 1) When to use this

Typical workflow for **WGBS / Bismark CX_report**:

1) Split each replicate:
```bash
python 02_prepare_cx_context_split.py -i rep1.CX_report.txt -p rep1 -d meth_data
python 02_prepare_cx_context_split.py -i rep2.CX_report.txt -p rep2 -d meth_data
python 02_prepare_cx_context_split.py -i rep3.CX_report.txt -p rep3 -d meth_data
```

2) Merge replicates per context:
```bash
python 03_prepare_cx_replicate_merge.py -o meth_data/merged_CG.CX  meth_data/rep1_CG.CX  meth_data/rep2_CG.CX  meth_data/rep3_CG.CX
python 03_prepare_cx_replicate_merge.py -o meth_data/merged_CHG.CX meth_data/rep1_CHG.CX meth_data/rep2_CHG.CX meth_data/rep3_CHG.CX
python 03_prepare_cx_replicate_merge.py -o meth_data/merged_CHH.CX meth_data/rep1_CHH.CX meth_data/rep2_CHH.CX meth_data/rep3_CHH.CX
```

---

## 2) Requirements

- Python ≥ 3.8
- No third-party dependencies (standard library only) fileciteturn56file0L55-L70

---

## 3) Input requirements (important)

### 3.1 Sorted inputs (mandatory)
Each input replicate file **must be sorted** by `(chrom, pos)` in ascending order. fileciteturn56file0L17-L23

Recommended sorting command (version-aware chromosome sort):
```bash
LC_ALL=C sort -k1,1V -k2,2n input.CX > input.sorted.CX
```

> Tip: if your chromosomes are like `Chr01, Chr02, ...`, `-V` works well.  
> If your chromosome naming is unusual, make sure the sort order is consistent across all replicates.

### 3.2 Column types
- `pos` must be an integer
- `col3` and `col4` must be integers fileciteturn56file0L84-L103  
The script casts them using `int(...)` and will error on non-integer values.

### 3.3 Missing positions across replicates
The output includes the **union** of positions across replicates:
- If a position exists in only one replicate, it is still written, and the sum uses only the present replicate(s). fileciteturn56file0L162-L192

---

## 4) Output

A 4-column file:

```
chrom   pos   sum_col3   sum_col4
```

Gzip output is enabled by using an output filename ending with `.gz`:

```bash
python 03_prepare_cx_replicate_merge.py -o merged_CG.CX.gz rep1_CG.CX.gz rep2_CG.CX.gz rep3_CG.CX.gz
```

---

## 5) Command-line arguments (checked against the script)

### Required
- `-o / --out FILE`  
  Output file path. If it ends with `.gz`, gzip output is produced. fileciteturn56file0L104-L115

- `inputs` (positional, **2 or 3** files)  
  Input replicate files (2 or 3). Each line must contain at least 4 columns:
  `chrom pos col3 col4`. Extra columns are ignored. fileciteturn56file0L10-L16 fileciteturn56file0L116-L121

> The script will exit if you provide fewer than 2 or more than 3 inputs. fileciteturn56file0L128-L132

---

## 6) Examples

### 6.1 Merge two replicates (plain text)
```bash
python 03_prepare_cx_replicate_merge.py -o merged_CG.CX rep1_CG.CX rep2_CG.CX
```

### 6.2 Merge three replicates (gzip input + gzip output)
```bash
python 03_prepare_cx_replicate_merge.py -o merged_CG.CX.gz rep1_CG.CX.gz rep2_CG.CX.gz rep3_CG.CX.gz
```

---

## 7) Troubleshooting

### “ERROR: please provide 2 or 3 input files.”
You passed the wrong number of replicate files. Provide exactly 2 or 3. fileciteturn56file0L128-L132

### “Bad line (need >=4 columns)”
One of your lines has fewer than 4 columns. Ensure the first 4 columns exist for every non-comment line. fileciteturn56file0L84-L95

### Output looks incomplete or inconsistent
Most commonly caused by **unsorted inputs**. Re-sort all inputs using the same sorting rule, then re-run the merge. fileciteturn56file0L17-L23

---

## 8) Help

```bash
python 03_prepare_cx_replicate_merge.py -h
```

(Argparse `prog` is `bts-cx-merge`; help output may show that name, but you can always run it via `03_prepare_cx_replicate_merge.py`.) fileciteturn56file0L104-L111
