# 12_plot_histone_vs_expr_heatmap.py — Track vs Expression Heatmap (EN)

Repo script: `12_plot_histone_vs_expr_heatmap.py`  
(Internal docstring/prog names may still show `omicscanvas_histone_vs_expr_heatmap_EN.py`; run by the repo filename above.)

This script plots **ChIP-seq / ATAC-seq (or any binned track) matrices** as heatmaps after **binning genes by RNA expression**.

It is designed to work **directly with the matrices produced by Step2** (`05_compute_bam_to_gene_matrices.py`).

---

## 1) Required inputs

### 1.1 Matrix directory (Step2 output)
`--matrix-dir` points to a directory containing per-track matrices for each region.

Recommended standardized suffixes (current OmicsCanvas convention):
- `*_tss_matrix.tsv`
- `*_gene_profile_matrix.tsv`
- `*_tes_matrix.tsv`

Examples (prefix = `B_`):
- `B_H3K4me3_tss_matrix.tsv`
- `B_H3K4me3_gene_profile_matrix.tsv`
- `B_H3K4me3_tes_matrix.tsv`

If you still have old filenames, override with:
- `--suffix-tss`, `--suffix-gene`, `--suffix-tes`

> Matrix bin columns: the script will try to **sort matrix columns numerically** (e.g. `0..199`).  
> If columns are not numeric-like strings, it keeps the original order.

### 1.2 Expression table (`--expr`)
Provide a TSV where **row index is gene ID**.

Example:
```text
GeneID	rep1	rep2	rep3
AT1G01010	5.2	4.9	5.1
AT1G01020	0	0	0
```

- `--expr-cols` selects which expression columns to average (0-based indices **after** the index column).
- `--expr-name` sets the internal column name used in logs/plots (default `EXP`).

---

## 2) Track definition syntax (`--tracks`)

`--tracks` is comma-separated. Each item supports:

1) `NAME`  
Display name = `NAME`; file id = `<file-prefix>NAME`

2) `NAME:FILEID`  
Display name = `NAME`; file id = `<file-prefix>FILEID`

3) Replicates: `NAME:REP1+REP2+REP3`  
Each `REP*` becomes `<file-prefix>REP*` if it does not already start with `--file-prefix`.  
Replicates are averaged element-wise before plotting.

Examples:
```bash
--tracks "ATAC,H3K4me1,H3K4me3" --file-prefix "B_"
--tracks "ATAC:ATAC_rep1+ATAC_rep2,H3K4me3" --file-prefix "B_"
```

---

## 3) Gene regions (`--gene-types`) and axis templates

`--gene-types` controls which region(s) are plotted:
- `TSS`  : around transcription start site
- `gene` : upstream + scaled gene body + downstream
- `TES`  : around transcription end site

Example:
```bash
--gene-types TSS,gene,TES
```

Axis ticks/labels and vertical guide lines are derived from:
- `--distance`
- `--bins-start`, `--bins-end`
- `--bins-gene-st`, `--bins-gene-body`, `--bins-gene-en`

**Important:** these `--bins-*` must match the actual matrix column count for each region, otherwise the script errors.

---

## 4) Gene ID matching (expression ↔ matrices)

Use:
- `--gene-id-mode exact` (default): IDs must match exactly
- `--gene-id-mode strip_dot`: remove isoform suffix after `.` (e.g., `AT1G01010.1 → AT1G01010`)

If `strip_dot` produces multiple isoforms per gene, `--collapse-isoform` controls how to merge:
- `mean` (default)
- `max`
- `first`

Optional filter:
- `--isoform-suffix .1` → keep only matrix rows ending with `.1` (applied before collapsing).

---

## 5) Expression binning

Genes are split into two parts:
- expression `<= --zero-threshold` → `--none-bins` groups
- expression `>  --zero-threshold` → `--exp-bins` groups

Each group produces **one heatmap row**.

Within each group:
- `--stat mean` (default) or `--stat median`

Sorting direction:
- default: low → high
- `--descending`: high → low

---

## 6) Color scaling (vmin/vmax) and colormap

### 6.1 Automatic scaling (recommended)
`--scale-mode`:
- `quantile` (default): use `--quantiles` (e.g., `0.01,0.99`)
- `diverging`: symmetric around 0 (good for Z-score matrices)
- `ratio`: non-negative (good for raw coverage)

Optional shortcut:
- `--range X`
  - diverging → `[-X, +X]`
  - others    → `[0, X]`

If `--cmap` is not set, the script chooses a default based on `--scale-mode`:
- ratio → `viridis`
- diverging/quantile → `RdBu_r`

### 6.2 Manual scaling
- `--vmin` and `--vmax` override automatic inference (best practice: set both together).

---

## 7) Outputs

For each (track × gene-type) combination, one heatmap file is written:

`<outdir>/<out-prefix>_<track>_<gene-type>_histone_vs_expr_heatmap.<ext>`

Where:
- `<ext>` is `pdf` or `png` controlled by `--out-format`
- `<track>` is the display name from `--tracks`, but it will be **sanitized** for filenames (non `[A-Za-z0-9._-]` → `_`).

Plot controls:
- `--dpi` (PNG and PDF rasterization)
- `--ytick-step` (row tick spacing)
- `--no-border` (hide outer border)
- `--title` (optional title)

---

## 8) Example commands (name unified)

### Example A: Standard naming, plot all regions
```bash
python 12_plot_histone_vs_expr_heatmap.py \
  --matrix-dir caculate_matrix \
  --tracks "ATAC,H3K4me1,H3K4me3" \
  --file-prefix "B_" \
  --gene-types TSS,gene,TES \
  --expr FPKM_all.txt \
  --expr-cols 0,1,2 \
  --none-bins 10 \
  --exp-bins 90 \
  --distance 2000 \
  --bins-start 100 --bins-end 100 \
  --bins-gene-st 50 --bins-gene-body 100 --bins-gene-en 50 \
  --scale-mode quantile --quantiles 0.01,0.99 \
  --out-format png \
  --outdir out_heatmap \
  --out-prefix B
```

### Example B: Old filenames (override suffix)
```bash
python 12_plot_histone_vs_expr_heatmap.py \
  --matrix-dir caculate_matrix \
  --tracks "H3K4me3" --file-prefix "B_" \
  --gene-types gene \
  --suffix-gene _gene_body_matrix.txt \
  --expr FPKM_all.txt \
  --out-prefix legacy
```

---

## 9) Full argument list (checked against the script)

Required:
- `--matrix-dir`
- `--tracks`
- `--expr`

Track IO:
- `--file-prefix` (default: empty)

Regions and suffix:
- `--gene-types` (default: `gene`)
- `--suffix-tss` (default: `_tss_matrix.tsv`)
- `--suffix-gene` (default: `_gene_profile_matrix.tsv`)
- `--suffix-tes` (default: `_tes_matrix.tsv`)

Axis template:
- `--distance` (default: 2000)
- `--bins-start` (default: 100)
- `--bins-end` (default: 100)
- `--bins-gene-st` (default: 50)
- `--bins-gene-body` (default: 100)
- `--bins-gene-en` (default: 50)

Gene ID matching:
- `--gene-id-mode {exact,strip_dot}` (default: exact)
- `--collapse-isoform {mean,max,first}` (default: mean)
- `--isoform-suffix` (default: None)

Expression:
- `--expr-cols` (default: `0,1,2`)
- `--expr-name` (default: `EXP`)

Binning:
- `--none-bins` (default: 10)
- `--exp-bins` (default: 90)
- `--zero-threshold` (default: 0)
- `--descending`

Row statistic:
- `--stat {mean,median}` (default: mean)

Color scaling:
- `--scale-mode {quantile,diverging,ratio}` (default: quantile)
- `--quantiles` (default: `0.01,0.99`)
- `--range` (default: None)
- `--vmin`, `--vmax`
- `--cmap` (default: auto by scale-mode)

Plot/output:
- `--ytick-step` (default: auto)
- `--no-border`
- `--title`
- `--dpi` (default: 300)
- `--out-format {pdf,png}` (default: pdf)
- `--out-prefix` (default: out)
- `--outdir` (default: `.`)

---

## 10) Help
```bash
python 12_plot_histone_vs_expr_heatmap.py -h
```
