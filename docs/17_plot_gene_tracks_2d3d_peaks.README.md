# 17_plot_gene_tracks_2d3d_peaks.py

OmicsCanvas — **Single-gene stacked tracks** (2D / pseudo-3D) with optional **MACS2 peak shading** and optional **methylation tracks**.

This script draws, for one target gene/transcript:

- **Gene model** from a GFF3 (UTR/CDS/intron backbone)
- **Coverage tracks** from one or more BAMs (ChIP/Input/RNA/ATAC/etc.)
- Optional **peak shading** from MACS2 `narrowPeak` / `broadPeak` / BED (per BAM track)
- Optional **methylation tracks** from prefix-based TSVs (`__CG/__CHG/__CHH`)

> Note: the internal argparse `prog` name is `omicscanvas_gene_tracks_2d3d.py`, but the actual file in this repo is `17_plot_gene_tracks_2d3d_peaks.py`. Run the script by its filename (examples below).

---

## 1) Requirements

### 1.1 Python packages

- `python>=3.8`
- `numpy`, `pandas`, `matplotlib`
- `pysam` (required for BAM coverage)

Conda example:

```bash
conda install -c conda-forge numpy pandas matplotlib
conda install -c bioconda pysam
```

---

## 2) Quick start

### 2.1 Minimal 2D plot (BAM only)

```bash
python 17_plot_gene_tracks_2d3d_peaks.py \
  --mode 2d \
  --gff3 genome.gff3 \
  --gene GeneX \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec  "H3K27me3.bam,H3K36me3.bam;ATAC.bam;RNA.bam" \
  --name-spec "H3K27me3,H3K36me3;ATAC;RNA" \
  --out GeneX_2d.pdf
```

### 2.2 2D plot + peak shading (per track)

```bash
python 17_plot_gene_tracks_2d3d_peaks.py \
  --mode 2d \
  --gff3 genome.gff3 \
  --gene GeneX \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec  "G_H3K27me3.bam,B_H3K27me3.bam,R_H3K27me3.bam;G_RNA.bam,B_RNA.bam,R_RNA.bam" \
  --name-spec "G,B,R;G,B,R" \
  --peak-dir peaks \
  --peak-spec "G_H3K27me3.narrowPeak,B_H3K27me3.narrowPeak,R_H3K27me3.narrowPeak;,,." \
  --color-map "G:#2ca02c,B:#1f77b4,R:#d62728" \
  --share-group-ylim \
  --fig-x 12 --fig-y 7 \
  --out GeneX_2d_peaks.pdf
```

### 2.3 3D plot + methylation (after BAM tracks)

```bash
python 17_plot_gene_tracks_2d3d_peaks.py \
  --mode 3d \
  --gff3 genome.gff3 \
  --gene GeneX \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec  "H3K27me3.bam,H3K36me3.bam;ATAC.bam;RNA.bam" \
  --name-spec "H3K27me3,H3K36me3;ATAC;RNA" \
  --meth-dir CX_gene \
  --meth-spec "SRRxxx__GeneX" \
  --meth-layout combined \
  --meth-position after \
  --fig-x 14 --fig-y 9 \
  --track-width 0.80 --track-height 0.12 \
  --x-off 0.10 --y-off 0.60 --layer-gap 0.05 \
  --out GeneX_3d_with_meth.png
```

---

## 3) Input specification

### 3.1 GFF3 gene model

Required:

- `--gff3`: GFF3 annotation file
- `--gene`: target **mRNA/transcript ID** to plot
- `--distance`: flank length around the gene (bp). The plotted window is:
  `[(mrna_start - distance), (mrna_end + distance)]`

GFF3 attribute keys (when extracting IDs):

- `--gene-keywords` (default: `ID`)  
  Attribute key used to read the mRNA/transcript ID from the `mRNA` feature.
- `--exon-keywords` (default: `Parent`)  
  Attribute key used to link exon/CDS/UTR features back to the transcript.

**Important**: this script is transcript-centric — it expects an `mRNA` record for the given ID. If your GFF3 uses different feature names/keys, change the keyword options accordingly.

### 3.2 BAM tracks (`--bam-spec`, `--name-spec`)

Required:

- `--bam-spec`: BAM files arranged by **groups**
- `--name-spec`: labels for BAM tracks (same structure as `--bam-spec`)
- `--bam-dir`: directory containing BAMs (default: `.`)

**Layer syntax (shared across multiple arguments)**

- `;` separates **groups** (each group becomes one row of tracks)
- `,` separates **tracks** inside a group

Example:

```txt
--bam-spec  "A.bam,B.bam;Input.bam" \
--name-spec "H3K27me3,H3K36me3;Input"
```

**Shape rule**: `--bam-spec` and `--name-spec` must have identical `;` group count and identical number of `,` items within each group.

**Do not use empty items** in `--bam-spec` / `--name-spec` (e.g. `A.bam,,C.bam`). Empty items are not preserved for these arguments and will break shape matching.

Coverage normalization:

- `--reads-length` (default: `150`) is used for scaling
- `--seq-type {paired,single}` (default: `paired`) affects the normalization factor

### 3.3 Optional peak shading (`--peak-spec`)

Peak shading overlays semi-transparent spans (`axvspan`) where peaks overlap the plotting window.

Optional:

- `--peak-spec`: peak files per BAM track (MACS2 `narrowPeak` / `broadPeak` / BED; **must match** the shape of `--name-spec`)
- `--peak-dir`: directory containing peak files (default: `.`)
- `--peak-placeholder`: placeholder token indicating “no peak file” (default: `.`)
- `--peak-color`: shading color (default: `black`)
- `--peak-alpha`: shading transparency (default: `0.18`)

How to specify missing peaks (per track):

- Explicit placeholder: `file1.narrowPeak,.,file3.narrowPeak`
- Empty item: `file1.narrowPeak,,file3.narrowPeak` (empty becomes placeholder automatically)

Accepted placeholders (treated as “no peak file”): empty string, `.`, `NA`, `N/A`, `NONE`, `null`, `-` (case-sensitive for some tokens; best practice is to use `.`).

**Compressed peaks**: peak files ending with `.gz` are supported.

### 3.4 Optional methylation tracks (`--meth-*`)

If enabled, the script loads methylation TSVs by **prefix**.

Inputs:

- `--meth-dir`: directory containing methylation TSVs
- `--meth-spec`: prefixes arranged using the same `;`/`,` layer syntax

For each prefix `<prefix>`, the script expects:

- `<prefix>__CG.tsv`
- `<prefix>__CHG.tsv`
- `<prefix>__CHH.tsv`

TSV format:

- Either **4 columns**: `name  po  me  al`
- Or **3 columns**: `po  me  al`
- Header is allowed (case-insensitive column names `name/po/me/al` are recognized)

Interpretation:

- `po` is used as an **index within the plotting window** `[0, window_len)`  
  (window start = `mrna_start - distance`)
- Plotted value = `ratio = me / al` (with `al=0` treated as 0)

Layout:

- `--meth-layout {combined,separate}` (default: `combined`)
  - `combined`: overlay CG/CHG/CHH in one axis
  - `separate`: draw three axes per prefix (CHH/CHG/CG)
- `--meth-position {before,after}` (default: `after`): position relative to BAM rows

Useful flags:

- `--meth-ylabel`: y-axis label for methylation (default: `mC`)
- `--meth-hide-prefix`: do not print prefix text inside the methylation axis
- `--meth-drop-me-zero`: drop rows where `me==0` **before** plotting (recommended for sparse data)

---

## 4) Styling and scaling

### 4.1 Colors (BAM tracks)

- `--cmap` (default: `Set2`): colormap for BAM tracks (auto-assign colors by label)

Two ways to enforce consistent colors:

**A) Label → color mapping (recommended)**  
`--color-map "LABEL:COLOR,LABEL:COLOR,..."` (also supports `=` and `;`)

Example:
```txt
--color-map "G:#2ca02c,B:#1f77b4,R:#d62728"
```

**B) Explicit per-track colors**  
`--color-spec` must match the `;`/`,` shape of `--name-spec`

Example:
```txt
--color-spec "#2ca02c,#1f77b4,#d62728;#2ca02c,#1f77b4,#d62728"
```

Color priority (highest → lowest):

1. `--color-spec`
2. `--color-map`
3. `--cmap`

### 4.2 Y-axis controls (BAM)

- `--ylim "ymin,ymax"`: global y-range for **all** BAM tracks
- `--group-ylims "ymin,ymax;ymin,ymax;..."`: one y-range per `;` group
- `--share-group-ylim`: auto-compute one shared y-range within each group

Rule of thumb:

- Use `--share-group-ylim` when you want fair within-group comparison.
- Use `--ylim` when you want identical scaling across all tracks.

### 4.3 Y-axis controls (methylation)

- `--meth-ylim "ymin,ymax"`: global y-range for all methylation tracks
- `--meth-group-ylims "ymin,ymax;ymin,ymax;..."`: per `;` group y-range
- `--meth-share-group-ylims`: auto-compute one shared ylim within each methylation group

---

## 5) 3D mode layout parameters (`--mode 3d`)

These arguments control the pseudo-3D stacking geometry (all are **figure-fraction** coordinates, 0..1):

- `--track-width`, `--track-height`: base panel size
- `--x-off`, `--y-off`: depth offsets (effective shift = factor × track width/height)
- `--base-left`, `--top-start`: anchor location
- `--layer-gap`: extra vertical gap between `;` groups

Tuning order:

1. Increase `--fig-x` / `--fig-y` first (most common fix)
2. Then adjust `--track-width` / `--track-height`
3. Finally tune `--x-off` / `--y-off` and `--layer-gap`

---

## 6) Troubleshooting

- **`pysam` not found / BAM plotting fails**  
  Install: `conda install -c bioconda pysam`

- **Empty tracks**  
  Ensure each BAM has an index (`.bai`) and chromosome naming matches between BAM and GFF3 (`chr1` vs `1`).

- **Gene not found**  
  The value passed to `--gene` must match the **mRNA/transcript ID** extracted from the GFF3 using `--gene-keywords` / `--exon-keywords`.

- **`--peak-spec` shape mismatch**  
  `--peak-spec` must match `--name-spec` exactly in `;` groups and `,` items. Use `.` or an empty item for tracks without peaks.

- **3D overlap / clipping**  
  Increase `--fig-x/--fig-y`, then shrink `--track-width/--track-height`, or reduce `--x-off/--y-off`, or increase `--layer-gap`.

---

## 7) Full CLI help

Run:

```bash
python 17_plot_gene_tracks_2d3d_peaks.py -h
```

