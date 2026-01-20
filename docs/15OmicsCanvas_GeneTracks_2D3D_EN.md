# OmicsCanvas — Single-gene Track Plot (2D/3D)

This script draws **stacked linear tracks** for a single gene: gene model (exons/introns) + multiple **BAM-derived coverage tracks** + optional **methylation tracks**. It supports **2D stacked** and **3D layered** rendering.

- Script (EN): `omicscanvas_gene_tracks_2d3d_EN.py`
- Script (CN): `omicscanvas_gene_tracks_2d3d_CN.py`

---

## 1) Inputs

### A) Genome annotation (GFF3)
Provide a GFF3 file for drawing the gene structure and for locating the target gene.

- `--gff3`: path to `*.gff3` (or `*.gff.gz` if your Python can read it)
- `--gene`: gene ID to plot (e.g., `Potri.001G055900.5.v3.0`)

If your GFF3 uses non-standard feature names, use:
- `--gene-keywords`: keywords/regex used to detect gene/transcript features (default: `gene|mRNA|transcript`)
- `--exon-keywords`: keywords/regex used to detect exon/CDS features (default: `exon|CDS`)

### B) BAM coverage tracks
Coverage is computed by reading BAM files around the gene region.

- `--bam-dir`: directory containing BAMs (default: `./bam`)
- `--bam-spec`: **track grouping spec** (see below)
- `--name-spec`: labels for each track in `--bam-spec` (same shape)

**Track grouping spec format**
- Use `;` to separate **groups** (rows / layers)
- Use `,` to separate **tracks** inside a group

Example:

```txt
--bam-spec  "H3K27me3.bam,H3K36me3.bam;ATAC.bam;RNA.bam" \
--name-spec "H3K27me3,H3K36me3;ATAC;RNA"
```

### C) Optional methylation tracks (prefix-based)
If provided, the script will load methylation TSVs by **prefix**.

- `--meth-dir`: folder containing methylation TSVs
- `--meth-spec`: same grouping syntax as BAM, but elements are **prefixes** (not full filenames)

The script expects per-prefix files in the form:
- `<prefix>__CG.tsv`
- `<prefix>__CHG.tsv`
- `<prefix>__CHH.tsv`

Example:

```txt
--meth-dir CX_gene \
--meth-spec "SRR9321764_CHH"  
```

> If your methylation TSV uses a different naming convention, change `load_meth_prefix()` inside the script.

---

## 2) Output

- `--out`: output figure path (recommended: `.png` for quick viewing or `.pdf` for vector)

Notes:
- **Large multi-track figures** can be huge as PDF. Prefer `png` for very large layouts.
- Use `--fig-x/--fig-y` to control the canvas size.

---

## 3) Key parameters

### Rendering
- `--mode {2d,3d}`: choose 2D stacked or 3D layered view
- `--fig-x`, `--fig-y`: figure size in inches
- `--cmap`: colormap for track heatmaps (default: `RdBu_r`)

### Genomic window / bins
- `--distance`: flank length around the gene (bp). For example `2000` means `-2000bp .. +2000bp` around TSS/TES.

### Coverage normalization
- `--reads-length`: approximate fragment length for extension (default: 200)
- `--seq-type {single,paired}`: how to interpret reads (default: `paired`)

### Y-axis controls (2D/3D)
- `--ylim "min,max"`: global y-limit for all tracks (e.g. `0,50`)
- `--group-ylims`: per-group y-limits using `;` (e.g. `0,10;0,50;0,5`)
- `--share-group-ylim`: share one y-limit per group

### Methylation options
- `--meth-layout {combined,separate}`: combined = CG/CHG/CHH merged; separate = 3 subtracks
- `--meth-ylabel`: y-axis label for methylation (default: `Methylation`)
- `--meth-hide-prefix`: hide prefix text in methylation labels
- `--meth-drop-me-zero`: drop sites with `me==0` before plotting (recommended)

### 3D layout (only used in `--mode 3d`)
These control the 3D “panel” placement.
- `--base-left`, `--top-start`
- `--track-width`, `--track-height`
- `--x-off`, `--y-off`
- `--layer-gap`

If you don’t want to tune these, start with defaults and only adjust figure size.

---

## 4) Examples

### Example 1 — 2D stacked track plot
```bash
python omicscanvas_gene_tracks_2d3d_EN.py \
  --mode 2d \
  --gff3 genome.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "H3K27me3.bam,H3K36me3.bam;ATAC.bam;RNA.bam" \
  --name-spec "H3K27me3,H3K36me3;ATAC;RNA" \
  --meth-dir CX_gene \
  --meth-spec "SRR9321764" \
  --meth-layout combined \
  --fig-x 12 --fig-y 6 \
  --out gene_tracks_2d.png
```

### Example 2 — 3D layered view
```bash
python omicscanvas_gene_tracks_2d3d_EN.py \
  --mode 3d \
  --gff3 genome.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "H3K27me3,H3K36me3;ATAC;RNA" \
  --name-spec "H3K27me3,H3K36me3;ATAC;RNA" \
  --fig-x 14 --fig-y 7 \
  --out gene_tracks_3d.png
```

---

## 5) Troubleshooting

- **No reads / empty tracks**: check BAM index (`.bai`) exists and region names match the GFF3 (chr naming).
- **Gene not found**: adjust `--gene-keywords/--exon-keywords`, or pass the transcript ID if your file is transcript-centric.
- **Layout overlaps**: increase `--fig-x/--fig-y`, then tune `--track-width/--track-height` (3D) or `--ylim`.
