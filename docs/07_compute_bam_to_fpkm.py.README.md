# 07_compute_bam_to_fpkm.py

OmicsCanvas utility — compute **strict FPKM / TPM** from a **coordinate-sorted + indexed BAM** using `pysam`.

This script is meant to work *after* Step 1 (GFF3 → BED + length table):

- BED intervals (at least 4 columns): `chrom  start  end  gene_id`
- Length table (2 columns, no header): `gene_id<TAB>length_bp` fileciteturn64file0L61-L68 fileciteturn64file7L80-L91

It then counts reads (SE) or fragments (PE) per interval under strict filters, and outputs:

`FPKM, TPM, counts, length_bp` fileciteturn64file0L65-L68 fileciteturn64file8L39-L52

---

## 1) What “strict” means in this script

### 1.1 Read-level filters (default strict)
A read is kept only if it is: mapped, not secondary, not supplementary, not duplicate (unless enabled), not QC-fail (unless enabled), and `MAPQ >= --min-mapq`. fileciteturn64file6L14-L35

Options:
- `--min-mapq INT` (default 0) fileciteturn64file0L72-L79
- `--include-duplicates` (default: exclude duplicates) fileciteturn64file0L75-L76
- `--include-qcfail` (default: exclude QC-fail) fileciteturn64file0L77-L78

### 1.2 Paired-end mode counts FRAGMENTS
In PE mode:
- **Library size** is counted as the number of fragments by counting only `read1` (one per template). fileciteturn64file6L110-L121
- **Interval counts** are fragments deduplicated by `query_name` within each interval (so a fragment is counted once even if both mates overlap). fileciteturn64file0L14-L34

Proper-pair strictness:
- By default, PE mode requires **proper pairs** and both mates mapped. fileciteturn64file6L38-L51
- To allow improper pairs, use `--allow-improper-pair` (this disables the proper-pair requirement). fileciteturn64file4L1-L4

> Note: if intervals overlap, the same fragment can be counted for multiple intervals (this is standard “overlap counting” behavior).

---

## 2) Coordinate conventions

- BED is expected to be **0-based, end-exclusive** (standard BED). fileciteturn64file0L55-L58
- If your “BED” is actually **1-based inclusive**, enable:
  - `--bed-one-based`  
  The script converts `[st, en]` → `[st-1, en)` (end unchanged). fileciteturn64file6L57-L77 fileciteturn64file4L6-L7

---

## 3) Output

`--out` writes a TSV with:

- index: `gene_id`
- columns: `FPKM`, `TPM`, `counts`, `length_bp` fileciteturn64file0L67-L68 fileciteturn64file8L50-L52

Formulas:
- `FPKM = 1e9 * counts / (length_bp * N)` where `N` is strict library size after filtering fileciteturn64file2L20-L40
- `TPM` is computed from `RPK = counts / (length_bp/1000)` and scaled to 1e6 fileciteturn64file8L42-L47

---

## 4) Requirements

- Python ≥ 3.8
- `pandas`
- `pysam` (required at runtime; imported lazily so `-h` still works without it) fileciteturn64file4L37-L43

Install:
```bash
pip install pandas pysam
```

---

## 5) Usage

### 5.1 Single-end (SE)
```bash
python 07_compute_bam_to_fpkm.py \
  -b genes.bed \
  -l gene_len.tsv \
  -m sample.bam \
  -o sample.fpkm_tpm.tsv \
  --mode se \
  --min-mapq 10
```

### 5.2 Paired-end (PE, strict proper pairs; default behavior)
```bash
python 07_compute_bam_to_fpkm.py \
  -b genes.bed \
  -l gene_len.tsv \
  -m sample.bam \
  -o sample.fpkm_tpm.tsv \
  --mode pe \
  --min-mapq 10
```

### 5.3 Paired-end but allow improper pairs
```bash
python 07_compute_bam_to_fpkm.py \
  -b genes.bed \
  -l gene_len.tsv \
  -m sample.bam \
  -o sample.fpkm_tpm.tsv \
  --mode pe \
  --allow-improper-pair
```

### 5.4 Auto-detect SE/PE (default)
```bash
python 07_compute_bam_to_fpkm.py \
  -b genes.bed \
  -l gene_len.tsv \
  -m sample.bam \
  -o sample.fpkm_tpm.tsv
```
Mode auto-detection inspects up to 200 mapped reads and chooses `pe` if any is paired. fileciteturn64file4L17-L30

---

## 6) Command-line arguments (complete)

Required:
- `-b/--bed` : BED with ≥4 columns (chrom, start, end, gene_id) fileciteturn64file0L61-L63
- `-l/--length` : length table (gene_id, length_bp) fileciteturn64file0L63-L64
- `-m/--bam` : BAM (sorted + indexed) fileciteturn64file0L65-L66
- `-o/--out` : output TSV fileciteturn64file0L67-L68

Mode:
- `--mode {auto,se,pe}` (default: `auto`) fileciteturn64file0L70-L71

Filtering:
- `--min-mapq INT` (default 0) fileciteturn64file0L72-L74
- `--include-duplicates` fileciteturn64file0L75-L76
- `--include-qcfail` fileciteturn64file0L77-L78

Paired-end strictness:
- `--require-proper-pair` (present in CLI; effectively enabled by default) fileciteturn64file4L1-L2
- `--allow-improper-pair` (disables the proper-pair requirement) fileciteturn64file4L3-L4

Coordinate conversion:
- `--bed-one-based` : convert 1-based inclusive BED to standard BED fileciteturn64file4L6-L7

Performance / logging:
- `--threads INT` (default 1) : threads for BAM BGZF decompression in pysam fileciteturn64file4L9-L10
- `--progress INT` (default 1000) : print progress every N intervals (0 disables) fileciteturn64file4L11-L12

---

## 7) Troubleshooting

### “BAM index not found or invalid”
The script requires a `.bai` and will raise an error if missing. fileciteturn64file4L48-L53  
Fix:
```bash
samtools index sample.bam
```

### “pysam is required”
Install:
```bash
pip install pysam
```
(raised by the script) fileciteturn64file4L39-L43

### “Library size N <= 0 after filtering”
Your filters removed all reads/fragments. Lower `--min-mapq`, or allow duplicates/QC-fail if appropriate. fileciteturn64file2L39-L41

### “Sum of RPK <= 0”
All counts are zero after filtering (or lengths are invalid). Check input files and filters. fileciteturn64file8L45-L47

---

## 8) Help
```bash
python 07_compute_bam_to_fpkm.py -h
```

(Help header may show `prog=omicscanvas_bam_to_fpkm.py`, but you can always run the repo script by filename.) fileciteturn64file0L42-L44
