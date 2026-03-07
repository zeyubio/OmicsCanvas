# 04_prepare_extract_gene_methylation.py

OmicsCanvas Step 4 — **Extract single-gene methylation (CX) records** for plotting.

This script **does not call methylation**. It only subsets your existing **CX coverage files**
(`ch  po  me  al`) by gene coordinates (± distance) and converts genomic coordinates into a
**strand-aware relative position** so that all genes align in the 5′→3′ direction. fileciteturn58file0L6-L25

It is typically used to generate the per-gene TSV inputs required by OmicsCanvas
single-gene plotting tools (gene track plots / circle plots), which expect files like:
`<prefix>__<gene>__CG.tsv` / `__CHG.tsv` / `__CHH.tsv`.

---

## 1) Requirements

- Python ≥ 3.8
- `numpy`, `pandas` fileciteturn58file0L31-L36

Install (pip):
```bash
pip install numpy pandas
```

---

## 2) Inputs

### 2.1 GFF3 annotation (required)

Arguments:
- `-g/--gff3` (required): GFF3 annotation file fileciteturn58file0L54-L63
- `--feature-type` (default: `mRNA`): which GFF3 feature type provides the target coordinates fileciteturn58file0L64-L71  
  Common values: `gene`, `mRNA`, `transcript`.
- `--attr-key` (default: `ID`): attribute key used to extract the identifier from the 9th column fileciteturn58file0L73-L80

ID extraction rule:
- The script searches the attribute field with regex: `<attr-key>=VALUE` and uses `VALUE` as the feature ID. fileciteturn58file0L73-L80

### 2.2 CX coverage files (required)

Arguments:
- `--cx-dir` (required): directory containing CX files. fileciteturn58file0L83-L90
- `--cx-suffix` (default: `.CX`): CX filename suffix; can be `.CX.gz` if gzipped. fileciteturn58file0L83-L90
- `--contexts` (default: `CG,CHG,CHH`): comma-separated contexts. fileciteturn58file0L83-L90

**Naming rule (important):**
For each SRR (or sample ID) and each context, the script loads:

`<cx-dir>/<SRR>_<context><cx-suffix>` fileciteturn58file0L45-L50

Example (default suffix `.CX`):
- `meth_data/SRR9321764_CHH.CX` fileciteturn58file1L34-L41

**File format (TSV; header not required):**
```
ch    po    me    al
```
- `ch`: chromosome (string)
- `po`: genomic position (integer; commonly 1-based)
- `me`: methylated counts (float/int)
- `al`: total counts / depth (float/int) fileciteturn58file0L14-L20 fileciteturn58file1L42-L50

> Tip (pipeline compatibility):  
> If you use `02_prepare_cx_context_split.py`, it produces files named `prefix_CG.CX`, etc.  
> To use them here, make sure your files follow `<SRR>_<context>.CX` (rename or symlink), and set `--srr` to the same prefix.

---

## 3) Target gene selection

You must provide **either** `--gene` **or** `--gene-list`. If both are provided, `--gene-list` takes precedence. fileciteturn58file0L146-L156

- `--gene`: single target ID (exact match to the extracted GFF3 ID)
- `--gene-list`: file with one ID per line (`#` comment lines are ignored); duplicates are removed while keeping order. fileciteturn58file0L157-L176

---

## 4) Relative coordinate convention (critical)

Let the extracted window be:
- `po1 = start - distance`
- `po2 = end + distance` fileciteturn58file1L55-L63

Relative coordinate is computed as:
- **‘+’ strand:** `rel_po = po - po1`
- **‘-’ strand:** `rel_po = po2 - po` fileciteturn58file0L22-L25 fileciteturn58file3L13-L18

So `rel_po` always increases along the transcription direction (5′→3′).

Window clipping:
- `po1` is clipped to ≥ 1 (so start-distance never becomes 0/negative). fileciteturn58file0L190-L201

---

## 5) Outputs

For each `(prefix, gene, context)` the script writes one TSV:

`<outdir>/<prefix>__<gene>__<context><out-suffix>` fileciteturn58file1L69-L78 fileciteturn58file3L85-L89

Default:
- `--outdir` = `single_gene`
- `--out-suffix` = `.tsv` fileciteturn58file0L118-L138

Columns (with header):
- `name` (gene ID)
- `po` (relative coordinate)
- `me` (methylated)
- `al` (total) fileciteturn58file0L27-L29

Filename safety:
- `prefix` and `gene` are sanitized so that any non `[A-Za-z0-9._-]` characters become `_`. fileciteturn58file0L178-L184

Overwrite behavior:
- If the output file exists and **`--overwrite` is NOT set**, the script will **skip** writing that file. fileciteturn58file3L88-L89

---

## 6) Command-line arguments (complete; checked against the script)

### Required
- `-g/--gff3 FILE` fileciteturn58file0L54-L63
- `--cx-dir DIR` fileciteturn58file0L83-L90
- `--srr LIST` : comma-separated SRR/sample IDs fileciteturn58file0L99-L106  
- One of: `--gene ID` or `--gene-list FILE` fileciteturn58file0L109-L126

### Optional (annotation)
- `--feature-type STR` (default: `mRNA`) fileciteturn58file0L64-L71
- `--attr-key STR` (default: `ID`) fileciteturn58file0L73-L80

### Optional (CX files)
- `--cx-suffix STR` (default: `.CX`) fileciteturn58file0L83-L90
- `--contexts STR` (default: `CG,CHG,CHH`) fileciteturn58file0L83-L90

### Optional (sample naming)
- `--prefix LIST` : optional comma-separated output prefixes
  - If not set: prefixes = SRR list
  - If shorter than SRR list: missing prefixes are auto-filled by the remaining SRR IDs fileciteturn58file3L38-L43

### Optional (window & IO)
- `--distance INT` (default: `1000`) fileciteturn58file0L127-L138
- `-o/--outdir DIR` (default: `single_gene`) fileciteturn58file0L118-L124
- `--out-suffix STR` (default: `.tsv`) fileciteturn58file0L118-L126
- `--chunksize INT` (default: `2000000`)  
  Read CX files in chunks to reduce memory. Set `0` to read the whole file at once. fileciteturn58file0L133-L138
- `--overwrite` : overwrite output files if they already exist fileciteturn58file0L139-L144

---

## 7) Examples (name unified)

### 7.1 One gene, one sample, all contexts
```bash
python 04_prepare_extract_gene_methylation.py \
  --gff3 genome/annotation.gff3 \
  --cx-dir meth/meth_data \
  --srr SRR9321764 \
  --gene Potri.001G055900.5.v3.0 \
  --contexts CG,CHG,CHH \
  --distance 2000 \
  --outdir single_gene
```
(Example matches the old guide logic but uses the repo script name.) fileciteturn58file1L84-L95

### 7.2 Multiple samples + gene list + custom prefixes
```bash
python 04_prepare_extract_gene_methylation.py \
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
(Structure matches the old guide logic but uses the repo script name.) fileciteturn58file1L97-L109

### 7.3 Gzipped CX files
```bash
python 04_prepare_extract_gene_methylation.py \
  --gff3 genome/annotation.gff3 \
  --cx-dir meth/meth_data \
  --cx-suffix .CX.gz \
  --srr SRR9321764 \
  --gene YOUR_GENE \
  --contexts CHH
```
(Structure matches the old guide logic but uses the repo script name.) fileciteturn58file1L110-L117

---

## 8) Troubleshooting

### “No target IDs found in GFF3”
- Check `--feature-type` (gene vs mRNA)
- Check `--attr-key` matches your attribute field (e.g., `ID`, `gene_id`, `transcript_id`) fileciteturn58file0L64-L80

### Missing CX file warnings
The script prints warnings and skips missing CX files. Ensure exact naming:
`<SRR>_<context><cx-suffix>`. fileciteturn58file3L67-L70

### Output is empty
- The CX file may not have cytosines inside the window
- Check whether your CX `po` is 1-based vs 0-based (this script treats it as whatever you provide; window is in GFF3 coordinates) fileciteturn58file0L14-L20

---

## 9) Help
```bash
python 04_prepare_extract_gene_methylation.py -h
```

(Argparse `prog` is `omicscanvas_extract_gene_methylation.py`; help output may show that name, but you can always run by `04_prepare_extract_gene_methylation.py`.) fileciteturn58file0L39-L42
