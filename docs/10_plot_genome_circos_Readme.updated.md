# 10_plot_genome_circos.py

Genome-wide Circos-like plot (bin-based) for multi-omics tracks (methylation CX, ChIP/ATAC BAM, RNA BAM) using **matplotlib only**.

## Features

- **Circos-like genome-wide view** using fixed **binning** (`--bin-size`)
- **Innermost ideogram** (chromosome band) + **uniform gene ticks** (green) from `--gene-bed`
- **Methylation (CX)** tracks (default: `spokes`)
- **Other omics BAM** tracks (default: `fill`)
- **RNA BAM** tracks (default: `fill`)
- **Y ticks**: only `0` and `max` for each track, displayed **once** at a fixed angle (default: right side)
- **Legend**: upper-right
- Chromosome selection by prefix + keyword filtering (no regex required)

**Track (ring) order**

- From inner (closest to ideogram) to outer rings:
  - `--cx-tracks` (in the given order)
  - `--bam-tracks` (in the given order)
  - `--rna-tracks` (in the given order)
- If a category is omitted, that category is simply skipped.

---

## Installation

```bash
pip install numpy pandas matplotlib pysam
```

> BAM files must be **sorted + indexed** (`.bai`).

---

## Quick Start

### Example (Poplar)

```bash
python 10_plot_genome_circos.py \
  --cx-dir  meth/meth_data \
  --bam-dir bam \
  --rna-dir bam \
  --cx-tracks  "SRR9321764_CG.CX,SRR9321764_CHG.CX,SRR9321764_CHH.CX" \
  --cx-names   "leaf_CG,leaf_CHG,leaf_CHH" \
  --bam-tracks "SRR8742441.sorted.bam,SRR8742377.sorted.bam,SRR8742376.sorted.bam,SRR8742375.sorted.bam,SRR8742374.sorted.bam,SRR8742373.sorted.bam" \
  --bam-names  "ATAC,H3K4me1,H3K4me3,H3K56ac,H3K36me3,H3K27me3" \
  --rna-tracks "SRR8742314.sorted.bam" \
  --rna-names  "leaf_RNA" \
  --fai genome/Ptrichocarpa_210_v3.0.fa.fai \
  --bin-size 100000 \
  --chrom-prefix Chr --max-chroms 5 \
  --gene-bed genome/Ptr_bed.txt \
  --omics-style fill --rna-style fill --cx-style spokes \
  --out-prefix results/poplar_reset
```

### Output

- `results/poplar_reset.circos.pdf` (default)
- or `results/poplar_reset.circos.png` if `--fig-format png`

---

## Input Formats

> **Header lines**: `--gene-bed` and `--cx-tracks` are read with `header=None`. So **do not include a header row** unless it is commented out with `#`.

### 1) Chromosome sizes

Choose one:

- `--fai genome.fa.fai` (recommended)
- `--chrom-sizes chrom.sizes` (2 columns: `chrom  length`, whitespace or tab)

### 2) Gene BED (`--gene-bed`)

At least 3 columns:

```text
chrom   start   end   ...
```

The script uses `midpoint = (start + end) / 2` and draws **equal-length green ticks** on the inner ideogram line.

> If gene ticks look like a green solid ring, downsample using `--gene-tick-stride`.

### 3) Methylation CX (`--cx-tracks`)

Each line (tab-separated):

```text
chrom   pos   meth_count   depth
```

Binning per window:
- `meth_bin = sum(meth_count)`
- `depth_bin = sum(depth)`
- `ratio = meth_bin / depth_bin` (bins with depth `< --min-bin-depth` set to 0)

### 4) BAM (`--bam-tracks`, `--rna-tracks`)

Requirements:
- sorted BAM
- indexed `.bai`

Counting:
- default PE mode counts only `read1` (`--pe-count read1`)
- default normalization `rpm` (`--norm rpm`)

---

## Parameter Reference (Detailed)

### A) Directories

- `--cx-dir`  
  Base directory for CX files. Relative paths in `--cx-tracks` are joined with this.

> Paths may be absolute or relative. `s3://...` and `gs://...` URIs are passed through as-is.

- `--bam-dir`  
  Base directory for non-RNA BAM (ChIP/ATAC/etc).

- `--rna-dir`  
  Base directory for RNA BAM.

---

### B) Gene ticks (innermost green ticks)

- `--gene-bed`  
  Gene BED file. If not provided, no gene ticks are drawn.

- `--gene-tick-color` (default `#00AA00`)  
  Color for gene ticks.

- `--gene-tick-len` (default `0.035`)  
  Tick length (radius units).

- `--gene-tick-lw` (default `0.6`)  
  Tick line width.

- `--gene-tick-stride` (default `1`)  
  Plot every N-th gene tick.  
  Example: `--gene-tick-stride 5` reduces density by 5×.

---

### C) Methylation tracks (CX)

- `--cx-tracks`  
  Comma-separated CX files: `a.CX,b.CX,c.CX`  
  Order = ring order (earlier = inner).

- `--cx-names`  
  Comma-separated labels (must match count of `--cx-tracks`).

- `--cx-colors`  
  Comma-separated colors (must match count of `--cx-tracks`).  
  If omitted, a default palette is used.

---

### D) Other omics BAM tracks

- `--bam-tracks`  
  Comma-separated BAM files.

- `--bam-names`  
  Comma-separated labels (must match count of `--bam-tracks`).

- `--bam-colors`  
  Comma-separated colors (must match count of `--bam-tracks`).

---

### E) RNA BAM tracks

- `--rna-tracks`  
  Comma-separated RNA BAM files.

- `--rna-names`  
  Comma-separated labels (must match count of `--rna-tracks`).

- `--rna-colors`  
  Comma-separated colors (must match count of `--rna-tracks`).

---

### F) Binning / CX filtering

- `--bin-size` (default `100000`)  
  Window size in bp. Common choices: `10000`, `50000`, `100000`.

- `--min-bin-depth` (default `10`)  
  For CX ratio, bins with depth `< this` are set to 0.

- `--cx-chunksize` (default `2000000`)  
  Pandas chunk size (rows) for CX reading.

---

### G) Chromosome sizes input

- `--fai`  
  FASTA index `.fai`.

- `--chrom-sizes`  
  Two-column file: `chrom length` (whitespace/tab).

---

### H) Chromosome selection

- `--chroms`  
  Explicit chromosome list: `Chr01,Chr02,...` (highest priority)

- `--chrom-prefix` (default `Chr`)  
  Keep contigs starting with this prefix.

- `--chrom-prefix-ignore-case`  
  Ignore case for prefix matching.

- `--exclude-keywords`  
  Exclude contigs containing these keywords (scaffolds/organelles).

- `--max-chroms` (default `25`)  
  Plot top N chromosomes by length (after filtering).

- `--min-chrom-len` (default `1`)  
  Minimum chromosome length to include (bp).

---

### I) BAM counting / filters

- `--pe-count` (`read1` or `all`, default `read1`)  
  For paired-end BAM, count only read1 (avoids double-count).

- `--include-duplicates`  
  Include duplicate reads (default: excluded).

- `--mapq` (default `1`)  
  Minimum MAPQ.

- `--norm` (`rpm` or `none`, default `rpm`)  
  `rpm` normalizes to reads-per-million.

---

### J) Track styles

- `--cx-style` (`spokes` or `fill`, default `spokes`)  
  CX display style.

- `--omics-style` (`spokes` or `fill`, default `fill`)  
  Style for non-RNA BAM tracks.

- `--rna-style` (`spokes` or `fill`, default `fill`)  
  Style for RNA tracks.

- `--fill-alpha` (default `0.7`)  
  Alpha transparency for `fill` tracks.

---

### K) Scaling / clipping

- `--scale-method` (`max`, `p95`, `p99`, `none`, default `p99`)  
  Scale used to map values to ring height.

- `--clip-mode` (`scale` or `none`, default `scale`)  
  Controls how `clip_val` is computed:
  - `scale`: `clip_val = scale * clip_mult` (stronger clipping when `clip_mult` < 1)
  - `none`: `clip_val = scale` (no extra multiplier; values are still capped at `scale`)

- `--clip-mult` (default `1.0`)  
  Clip threshold multiplier: `clip_val = scale * clip_mult`.

---

### L) Chromosome spacing

- `--gap-frac` (default `0.06` in v2_1)  
  Gap size between chromosomes as a fraction of chromosome bin count.

- `--gap-bins` (default `None`)  
  Fixed gap in bins (overrides `gap-frac`).

- `--no-gap-after-last`  
  Disable the extra gap after the last chromosome (default: enabled).

---

### M) Figure / layout

- `--fig-size` (default `10.0`)  
  Figure size in inches.

- `--dpi` (default `300`)  
  DPI for PNG output.

- `--fig-format` (`pdf` or `png`, default `pdf`)  

- `--ideogram-radius` (default `1.0`)  
  Outer radius of ideogram band.

- `--ideogram-width` (default `0.08`)  
  Width of ideogram band.

- `--chr-label-radius` (default `0.70`)  
  Radius for chromosome name labels (placed in the inner blank area).

- `--base-radius` (default `1.10`)  
  Radius of the first data track.

- `--track-height` (default `0.20`)  
  Height per track.

- `--track-gap` (default `0.13`)  
  Gap between tracks.

- `--rotate-deg` (default `0.0`)  
  Rotation of the whole plot. `0` means first chromosome starts at right.

- `--ytick-angle-deg` (default `0.0`)  
  Angle for displaying y ticks (0/max).  
  - `0` = right side (recommended)
  - `90` = top

- `--outline-lw` (default `0.8`)  
  Ideogram outline width.

- `--seg-lw` (default `0.25`)  
  Line width for spokes.

---

### N) Legend

- `--no-legend`  
  Disable legend.

- `--legend-fontsize` (default `10`)  
  Legend font size.

---

### O) Output

- `--out-prefix` (required)  
  Output prefix. Figure saved as `<out-prefix>.circos.<pdf|png>`.

---

## Tips

### Gene ticks too dense

```bash
--gene-tick-stride 5
```

### Chromosomes not selected as expected

- Check contig names in `.fai`
- Adjust:
  - `--chrom-prefix`
  - `--exclude-keywords`
  - or specify exact list with `--chroms`

### BAM is slow

Use larger `--bin-size` (e.g. `100000`) and fewer chromosomes (`--max-chroms`).

---

## License / Citation

Add your project license and citation here if needed.
