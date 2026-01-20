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

