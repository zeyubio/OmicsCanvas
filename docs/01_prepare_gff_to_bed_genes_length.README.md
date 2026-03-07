# 01_prepare_gff_to_bed_genes_length.py

OmicsCanvas Step 1 — Convert **GFF3 → BED6** and compute **CDS/exon length**.

This script exports coordinates for `gene` or `mRNA/transcript` features as BED6, and summarizes **CDS** or **exon** length.
It also auto-resolves common ID mismatches (no `--level` needed):

- **If BED IDs already match** CDS/exon `Parent` IDs → write length table directly.
- **If you build BED from genes** (`--bed-feature gene`) but CDS/exon `Parent` points to transcripts → the script builds a transcript→gene map from `mRNA/transcript` records (`ID=<tx>`, `Parent=<gene>`) and aggregates transcript lengths to genes. fileciteturn51file3 fileciteturn51file4

---

## 1) Requirements

- Python ≥ 3.8
- `pandas`

Install (pip):
```bash
pip install pandas
```

---

## 2) Input

### 2.1 GFF3 format requirements

- Standard **GFF3** (9 columns, tab-delimited).
- Must include:
  - BED feature records: `gene` or `mRNA/transcript`
  - Length feature records: `CDS` or `exon`
- Attribute column must include:
  - ID key (default `ID=...`) for the BED feature
  - Parent key (default `Parent=...`) for CDS/exon (and for transcript→gene mapping when needed)

**Multi-parent supported:** `Parent=a,b,c` will be split and exploded automatically. fileciteturn51file2L1-L8

---

## 3) Outputs

### 3.1 BED6 (`--out-bed`, default: `gene.bed`)

Format (no header):
```
chrom   start   end   ID   score   strand
```

Notes:
- The output filename default is `gene.bed`, **but the content depends on `--bed-feature`**.
  - Default `--bed-feature` is `mRNA` fileciteturn51file1L22-L26
  - So if you do not set `--bed-feature`, `gene.bed` will actually contain transcript/mRNA records.

### 3.2 Length table (`--out-length`, default: `feature_length.tsv`)

Two columns (no header):
```
ID<TAB>length
```

- Length is computed from CDS/exon segments using **1-based inclusive** coordinates.
- With `--merge-overlaps`, overlapping segments are merged before summing (recommended when overlaps exist). fileciteturn51file1L44-L45

---

## 4) Coordinate conventions

GFF3 coordinates are **1-based inclusive**. BED is **0-based half-open**.

- Default: the script writes coordinates as-is from GFF3.
- If you set `--bed-zero-based`, the script converts:
  - `start = start - 1`
  - `end` unchanged  
  This matches BED’s 0-based half-open convention. fileciteturn51file1L34-L39 fileciteturn51file3L6-L10

---

## 5) Command-line arguments (checked against the script)

### Required
- `-i/--gff3` : input GFF3 file path. fileciteturn51file1L14-L20

### Output
- `-o/--out-bed` (default: `gene.bed`) : output BED6 path. fileciteturn51file1L14-L20
- `-l/--out-length` (default: `feature_length.tsv`) : output length table path. fileciteturn51file1L17-L20

### Feature selection
- `--bed-feature {gene,mRNA,transcript}` (default: `mRNA`)  
  Feature type used to build BED6. fileciteturn51file1L22-L26
- `--length-feature {CDS,exon}` (default: `CDS`)  
  Feature type used to compute length. fileciteturn51file1L24-L26

### Attribute keys
- `--id-key` (default: `ID`)  
  Attribute key for feature ID (e.g. `ID=xxx`). fileciteturn51file1L28-L33
- `--parent-key` (default: `Parent`)  
  Attribute key for parent ID (e.g. `Parent=xxx`). fileciteturn51file1L30-L33

### Coordinate conversion
- `--bed-zero-based`  
  Convert GFF3 start to BED 0-based half-open (start-1, end unchanged). fileciteturn51file1L34-L39

### Other options
- `--score` (default: `10`)  
  Constant value for BED score column (5th column). fileciteturn51file1L41-L43
- `--merge-overlaps`  
  Merge overlapping CDS/exon segments before summing length (more accurate). fileciteturn51file1L44-L45
- `--min-match` (default: `10000`)  
  Minimum number of shared IDs between BED IDs and length IDs. If the intersection is smaller, the script aborts to prevent wrong mapping. Lower this for small datasets or when working with a subset. fileciteturn51file1L46-L51

---

## 6) Examples (name unified to this script)

### 6.1 BED from transcripts (mRNA), CDS length by transcript
```bash
python 01_prepare_gff_to_bed_genes_length.py   -i annotation.gff3   --bed-feature mRNA   --length-feature CDS   --bed-zero-based   -o transcript.bed   -l transcript_cds_length.tsv   --merge-overlaps   --min-match 1000
```

### 6.2 BED from genes, aggregate transcript CDS length to genes (auto mapping)
```bash
python 01_prepare_gff_to_bed_genes_length.py   -i annotation.gff3   --bed-feature gene   --length-feature CDS   --bed-zero-based   -o gene.bed   -l gene_cds_length.tsv   --merge-overlaps   --min-match 1000
```

### 6.3 Non-standard GFF3 attribute keys
If your GFF3 uses different keys (e.g. `gene_id=` instead of `ID=`):
```bash
python 01_prepare_gff_to_bed_genes_length.py   -i annotation.gff3   --bed-feature gene   --id-key gene_id   --parent-key Parent   --length-feature exon   --bed-zero-based   -o gene.bed   -l gene_exon_length.tsv   --min-match 100
```

---

## 7) Troubleshooting

### “ID mismatch … intersection < --min-match”
This means BED IDs and CDS/exon Parent IDs do not match at the level you selected.

Fixes:
- If CDS/exon `Parent` points to transcript IDs (common), and you want gene-level lengths:
  - Use `--bed-feature gene` (the script will auto aggregate tx→gene). fileciteturn51file4L1-L36
- If your GFF3 uses different attribute keys:
  - Set `--id-key` and/or `--parent-key`. fileciteturn51file4L41-L45
- If you run on a subset of genes:
  - Lower `--min-match` (e.g. 100 / 10).

### BED coordinates look “shifted”
- If you need standard BED (0-based half-open), add `--bed-zero-based`. fileciteturn51file1L34-L39

---

## 8) CLI help

```bash
python 01_prepare_gff_to_bed_genes_length.py -h
```
