# OmicsCanvas – Extract Gene Methylation (CX) Records

This document describes how to use **`omicscanvas_extract_gene_methylation_EN.py`** to extract
per-base methylation (CX) records for **one gene** or a **gene list**, producing a TSV that can be
used by downstream OmicsCanvas plotting scripts.

---

## 1. What this script does

- Reads **gene/transcript coordinates** from a **GFF3** file.
- Reads **CX methylation files** (precomputed) in the format: `ch  po  me  al`.
- For each target gene, extracts records in `[start-distance, end+distance]`.
- Converts genomic positions to **strand-aware relative coordinates** (so all genes align 5'→3').

It **does not** compute methylation; it only subsets your existing CX files.

---

## 2. Input requirements

### 2.1 GFF3
The script reads these columns:
- `seqid` (chromosome)
- `type` (feature type, e.g. mRNA or gene)
- `start`, `end`
- `strand` (+ or -)
- attributes (9th column)

You must specify:
- `--feature-type` (default: `mRNA`)
- `--attr-key` to extract the target ID from attributes (default: `ID`)

### 2.2 CX files
Expected naming under `--cx-dir`:

- `<SRR>_<context><cx-suffix>`

Default `cx-suffix` is `.CX`, so for example:
- `SRR9321764_CHH.CX`

Expected file format (no header by default):

```
ch\tpo\tme\tal
```

- `po` should be integer genomic coordinate (commonly 1-based)
- `me`, `al` can be float/int (counts)

---

## 3. Relative coordinate convention (important)

Let the extracted window be:
- `po1 = start - distance`
- `po2 = end + distance`

Relative coordinate is:
- `+` strand: `rel_po = po - po1`
- `-` strand: `rel_po = po2 - po`

So **rel_po always increases along transcription direction (5'→3')**.

---

## 4. Output

For each combination of (SRR/prefix, gene, context), the script writes:

- `<outdir>/<prefix>__<gene>__<context>.tsv`

Columns:

- `name`: gene ID
- `po`: relative coordinate
- `me`: methylated counts (numerator)
- `al`: total counts (denominator)

---

## 5. Common usage

### 5.1 Extract one gene from one SRR (all contexts)

```bash
python omicscanvas_extract_gene_methylation_EN.py \
  --gff3 genome/annotation.gff3 \
  --cx-dir meth/meth_data \
  --srr SRR9321764 \
  --gene Potri.001G055900.5.v3.0 \
  --contexts CG,CHG,CHH \
  --distance 2000 \
  --outdir single_gene
```

### 5.2 Extract a gene list from multiple SRRs

```bash
python omicscanvas_extract_gene_methylation_EN.py \
  --gff3 genome/annotation.gff3 \
  --cx-dir meth/meth_data \
  --srr SRR1,SRR2,SRR3 \
  --prefix WT,mut1,mut2 \
  --gene-list genes.txt \
  --contexts CHH \
  --distance 1000 \
  --outdir out_single_gene \
  --chunksize 2000000
```

### 5.3 If your CX files are gzipped

```bash
python omicscanvas_extract_gene_methylation_EN.py \
  --gff3 genome/annotation.gff3 \
  --cx-dir meth/meth_data \
  --cx-suffix .CX.gz \
  --srr SRR9321764 \
  --gene YOUR_GENE \
  --contexts CHH
```

---

## 6. Performance notes

- The script groups target genes by chromosome.
- For each chromosome, it reads only the **union range** covering all target genes on that chromosome.
- Use `--chunksize` to reduce memory; set `--chunksize 0` if you have enough RAM and want faster IO.

---

## 7. Troubleshooting

1) **No target genes found**
- Verify `--feature-type` (gene vs mRNA)
- Verify `--attr-key` matches your attribute field

2) **Missing CX file warnings**
- Ensure files are named exactly `<SRR>_<context><cx-suffix>`

3) Output is empty
- Your CX file may not contain records in that interval
- Check whether `po` is 1-based vs 0-based in your CX file

