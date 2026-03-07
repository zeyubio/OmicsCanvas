# 18_plot_gene_circle_plot_peaks.py

Single-gene **circle / ring** plot for OmicsCanvas with:
- gene structure ring (exon/CDS/UTR + strand arrow)
- optional **BAM coverage rings** (scaled coverage “spikes”)
- optional **MACS2 peak shading** aligned to BAM rings (narrowPeak/broadPeak/BED)
- optional **DNA methylation rings** (CG/CHG/CHH; combined or separate)

---

## 1. Quick start

### 1.1 Minimal (gene structure ring only)
```bash
python 18_plot_gene_circle_plot_peaks.py \
  --gff3 genome.gff3 \
  --gene GeneX \
  --out GeneX.circle.png
```

### 1.2 With BAM rings
```bash
python 18_plot_gene_circle_plot_peaks.py \
  --gff3 genome.gff3 \
  --gene GeneX \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec  "H3K27me3.bam,H3K36me3.bam;ATAC.bam" \
  --name-spec "H3K27me3,H3K36me3;ATAC" \
  --out GeneX.circle.pdf
```

### 1.3 With BAM rings + peak shading (recommended for enhancer context)
```bash
python 18_plot_gene_circle_plot_peaks.py \
  --gff3 genome.gff3 \
  --gene GeneX \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec  "H3K27ac.bam,ATAC.bam" \
  --name-spec "H3K27ac,ATAC" \
  --peak-dir peaks \
  --peak-spec "H3K27ac_peaks.narrowPeak,ATAC_peaks.narrowPeak" \
  --peak-alpha 0.18 \
  --out GeneX.circle_peaks.png
```

### 1.4 With methylation rings
```bash
python 18_plot_gene_circle_plot_peaks.py \
  --gff3 genome.gff3 \
  --gene GeneX \
  --distance 2000 \
  --meth-dir CX_gene \
  --meth-spec "Sample1__GeneX;Sample2__GeneX" \
  --circle-meth-layout combined \
  --out GeneX.circle_meth.png
```

---

## 2. Inputs and formats

### 2.1 GFF3 + target gene ID (required)
- `--gff3`: genome annotation in **GFF3** (tab-delimited, standard 9 columns)
- `--gene`: **mRNA / transcript ID** present in the GFF3
- `--distance`: flanking window in bp (default `1000`)

This tool is **transcript-centric**: it searches the `mRNA` feature (default type `mRNA`) and uses:
- `--gene-keywords` (default `ID`) to extract the mRNA ID from the GFF3 attributes
- `--exon-keywords` (default `Parent`) to link exon/CDS/UTR features back to the transcript

If your GFF3 uses different attribute keys, override these two options.

### 2.2 Layer syntax for grouped inputs
The following arguments use the same “layer” syntax:
- `--bam-spec`, `--name-spec`
- `--peak-spec`
- `--meth-spec`

Rules:
- `;` separates **groups**
- `,` separates **tracks within a group**

Example (2 groups; first has 2 tracks, second has 1 track):
```txt
--bam-spec  "A.bam,B.bam;C.bam"
--name-spec "A,B;C"
```

### 2.3 BAM rings (optional)
If `--bam-spec` is not set, BAM rings are skipped.

**Inputs**
- `--bam-dir`: directory containing BAMs (default `.`)
- `--bam-spec`: BAM filenames using the layer syntax
- `--name-spec`: labels aligned to `--bam-spec` (same structure)

**Coverage normalization**
Coverage is derived from `pysam.AlignmentFile.count_coverage()` and normalized by mapped reads (and by read length).  
Use:
- `--reads-length` (default `150`) to match your effective fragment length
- `--seq-type {paired,single}` (default `paired`) to match your library type

### 2.4 Peak shading (optional; requires BAM rings)
Peak shading is drawn as semi-transparent sectors on top of the corresponding BAM ring.

**Inputs**
- `--peak-dir`: directory containing peak files (default `.`)
- `--peak-spec`: peak filenames aligned to the BAM rings (same `;` and `,` shape as `--name-spec`)
- `--peak-placeholder`: placeholder meaning “no peaks for this ring” (default `.`)

Supported formats:
- MACS2 `narrowPeak`, `broadPeak`
- any BED-like file with at least 3 columns: `chrom  start  end`
- `.gz` and `.zip` are supported (zip uses the first non-directory member)

**Alignment constraint**
- You may only use `--peak-spec` when `--bam-spec/--name-spec` is provided.
- Peak entries must align 1-to-1 with the BAM rings. Use `.` (or an empty token between commas) to keep alignment.

Examples:
```txt
# 3 BAM rings; the middle ring has no peaks
--peak-spec "H3K27ac_peaks.narrowPeak,,ATAC_peaks.narrowPeak"

# Equivalent explicit placeholder
--peak-spec "H3K27ac_peaks.narrowPeak,.,ATAC_peaks.narrowPeak"
```

### 2.5 Methylation rings (optional)
If `--meth-spec` is set, `--meth-dir` is required.

**Prefix-based loading**
For each prefix `<prefix>` in `--meth-spec`, the tool loads:
- `<prefix>__CG.tsv`
- `<prefix>__CHG.tsv`
- `<prefix>__CHH.tsv`

**TSV formats**
Either with header (case-insensitive) or without header:
- 4 columns: `name  po  me  al`
- 3 columns: `po  me  al`

Interpretation:
- `po` is treated as the **0-based position within the plotting window** `[0, window_len)`
- plotted value is `ratio = me / al` (with `al==0` handled as 0)

If a `name` column exists, rows are filtered to `name == --gene` (with a warning if nothing matches).

**Layout**
- `--circle-meth-layout combined` (default): CG/CHG/CHH overlaid in one ring
- `--circle-meth-layout separate`: one ring per context

---

## 3. Command-line arguments

### 3.1 Required
- `--out` : output figure path (`.pdf`, `.png`, `.svg`)
- `--gff3`: GFF3 annotation file
- `--gene`: target transcript (mRNA ID) in GFF3

### 3.2 Gene model
- `--distance` (int, default `1000`): flank size in bp
- `--gene-keywords` (str, default `ID`): attribute key for mRNA ID
- `--exon-keywords` (str, default `Parent`): attribute key for exon/CDS parent

### 3.3 BAM rings (optional)
- `--bam-spec` (str): BAM layer spec (`;` groups, `,` tracks). If omitted, BAM rings are skipped.
- `--name-spec` (str): labels aligned to `--bam-spec`. Required when `--bam-spec` is provided.
- `--bam-dir` (str, default `.`): directory of BAM files
- `--reads-length` (int, default `150`): read/fragment length scaling factor for normalization
- `--seq-type {paired,single}` (default `paired`): affects normalization
- `--circle-bam-step` (int, default `10`): down-sample step for BAM spikes (larger = faster / sparser)
- `--bam-group-ylims` (default `auto`): per-group y-limits for BAM rings
  - `auto` (default): each ring scales independently
  - `MAX1;MAX2;...`: interpreted as `(0,MAX)` for each group
  - `MIN,MAX;MIN,MAX;...`: explicit limits
  - must match the number of `;` groups in `--bam-spec/--name-spec`
  - example: `--bam-group-ylims "50;50;200"`

### 3.4 Peak shading (optional; requires BAM rings)
- `--peak-spec` (str): peak files aligned to BAM rings (same shape as `--name-spec`)
- `--peak-dir` (str, default `.`): directory containing peak files
- `--peak-placeholder` (str, default `.`): placeholder meaning “no peak”
- `--peak-color` (str, default `black`): shading color
- `--peak-alpha` (float, default `0.18`): shading transparency
- `--peak-arc-points` (int, default `30`): polygon smoothness (larger = smoother, slower)

### 3.5 Methylation rings (optional)
- `--meth-dir` (str): directory of methylation TSVs (required if `--meth-spec` is provided)
- `--meth-spec` (str): prefix layer spec (`;` groups, `,` tracks)
- `--circle-meth-layout {combined,separate}` (default `combined`)
- `--meth-drop-me-zero` (flag): drop rows with `me==0` when loading methylation TSVs (often reduces noise)

### 3.6 Figure & ring geometry
- `--fig-x` (float, default `10.0`): figure width (inches)
- `--fig-y` (float, default `10.0`): figure height (inches)
- `--circle-ring-start` (float, default `1.10`): radius of the first ring outside the gene ring
- `--circle-ring-width` (float, default `0.20`): ring width (also the spike-length scale)
- `--circle-ring-gap` (float, default `0.08`): radial gap between successive rings

### 3.7 Frame & limits
- `--circle-frame` / `--circle-no-frame`: enable/disable the outer frame box (default: enabled)
- `--circle-frame-lw` (float, default `1.0`): frame line width
- `--circle-margin` (float, default `0.15`): extra margin for **auto** x/y limits
- `--circle-lim` (float, default `None`): manual half-range for x/y limits  
  Example: `--circle-lim 2.2` sets `xlim = ylim = [-2.2, 2.2]`

> Note: `--circle-lim` controls **plot limits only** (to prevent clipping), not per-ring scaling.
> For BAM ring y-scaling, use `--bam-group-ylims`.

---

## 4. Practical tips

- If the plot looks “too dense”: increase `--circle-bam-step` and/or increase `--fig-x/--fig-y`.
- If rings are clipped: increase `--circle-margin` or set `--circle-lim` explicitly.
- If you want fair comparison inside a sample: put tracks in the same `;` group and set `--bam-group-ylims`.
- For enhancer candidates: typical combination is H3K27ac + ATAC + peaks on both rings.

---

## 5. Troubleshooting

### 5.1 `pysam is required ...`
You enabled BAM rings (`--bam-spec`) but pysam is not installed.

Install:
```bash
conda install -c bioconda pysam
```

### 5.2 `--peak-spec requires --bam-spec/--name-spec`
Peak shading is only defined relative to BAM rings. Provide `--bam-spec` and `--name-spec`, or remove `--peak-spec`.

### 5.3 `Missing methylation file: <prefix>__CG.tsv`
Check that:
- `--meth-dir` is correct
- the filenames match `<prefix>__CG.tsv`, `<prefix>__CHG.tsv`, `<prefix>__CHH.tsv`

### 5.4 `No rows matched gene_id=...` warning for methylation
Your methylation TSVs include a `name` column, but its values do not match `--gene`.  
Either:
- use the transcript ID that matches `name`, or
- export methylation TSVs for the same gene ID you pass via `--gene`.

---

## 6. Reproducibility notes
- BAM coverage uses raw count_coverage + mapped-read normalization; results depend on BAM filtering and duplicates.
- Peaks are treated as simple intervals (chrom/start/end) and are intersected with the gene window `[mrna_start-distance, mrna_end+distance]`.
