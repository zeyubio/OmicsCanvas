# OmicsCanvas — Single-gene Circle Plot

This script renders a **circular** view of a single gene, with optional **BAM coverage rings** and **methylation rings**.

- Script (EN): `omicscanvas_gene_circle_plot_EN.py`
- Script (CN): `omicscanvas_gene_circle_plot_CN.py`

---

## 1) Inputs

### A) GFF3 + Gene ID
- `--gff3`: genome annotation (GFF3)
- `--gene`: gene/transcript ID
- `--distance`: flank length (bp)

### B) BAM rings
- `--bam-dir`: directory containing BAMs
- `--bam-spec`: `;` groups and `,` tracks (same syntax as the 2D/3D script)
- `--name-spec`: labels aligned to `--bam-spec`

### C) Methylation rings (optional)
- `--meth-dir`: directory of methylation TSVs
- `--meth-spec`: prefix list (same grouping syntax)

Expected prefix files:
- `<prefix>__CG.tsv`, `<prefix>__CHG.tsv`, `<prefix>__CHH.tsv`

---

## 2) Output
- `--out`: output path (recommend `.png` for convenience)

---

## 3) Key parameters

### Figure
- `--fig-x`, `--fig-y`: figure size in inches

### BAM processing
- `--reads-length`: fragment length extension (default 200)
- `--seq-type {single,paired}`

### Circle layout
- `--circle-bam-step`: down-sampling step for BAM points on the ring (smaller = denser)
- `--circle-ring-start`: inner radius to start drawing rings
- `--circle-ring-width`: width of each ring
- `--circle-ring-gap`: gap between rings

### Style
- `--circle-frame / --no-circle-frame`: draw outer frame
- `--circle-frame-lw`: frame line width
- `--circle-margin`: extra margin
- `--circle-lim`: manual color/height scaling for rings

### Methylation
- `--circle-meth-layout {combined,separate}`: combined = merge CG/CHG/CHH; separate = three rings
- `--meth-drop-me-zero`: drop `me==0` sites

---

## 4) Example
```bash
python omicscanvas_gene_circle_plot_EN.py \
  --gff3 genome.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "H3K27me3.bam,H3K36me3.bam;ATAC.bam" \
  --name-spec "H3K27me3,H3K36me3;ATAC" \
  --meth-dir CX_gene \
  --meth-spec "SRR9321764" \
  --circle-meth-layout combined \
  --fig-x 10 --fig-y 10 \
  --out gene_circle.png
```

---

## 5) Troubleshooting
- If the circle looks too dense: increase `--circle-bam-step` or increase figure size.
- If rings are clipped: increase `--circle-margin` or disable frame.
