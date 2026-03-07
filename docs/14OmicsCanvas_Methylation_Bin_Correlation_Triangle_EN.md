# OmicsCanvas — Methylation Bin Correlation (Triangle Heatmap)

Repo script: `14_methylation_bin_correlation_triangle.py`  
(Internal docstring name may still show `17_methylation_bin_correlation_triangle.py`; consider renaming for consistency.)

This tool:

1. **Bins** cytosine methylation “CX-like” files into fixed-size genomic bins (bp).
2. Computes per-bin methylation ratio: **sum(methylated) / sum(depth)**.
3. Intersects bins across samples (keeps only bins shared by all samples).
4. Computes **sample–sample correlation** (Pearson/Spearman/Kendall).
5. Plots a **triangular correlation heatmap** using **diamond cells** (one triangle only).

It supports **multi-panel** output via `|`-separated groups (e.g., `CG|CHG|CHH` stacked panels).

---

## 1) Input format

### 1.1 CX-like methylation file (required)

Tab-delimited, **no header** (header is allowed only if it is commented out with `#`).

Columns:

1. `chrom` (string)
2. `pos` (1-based integer)
3. `meth` (methylated count)
4. `depth` (coverage / depth)

Example:

```
Chr01   15072   1   1
Chr01   15073   0   1
Chr01   15347   8   8
Chr01   15348   4   4
```

Notes:
- Rows with `depth <= 0` are ignored.
- Bins are computed by: `bin_id = (pos-1) // bin_size`.

---

## 2) Outputs

For each panel `<group>`:

- `<out-prefix>.<group>.corr.tsv`  
  Correlation matrix (samples × samples), TSV.

Optional:
- `<out-prefix>.<group>.bin_matrix.tsv`  
  Per-bin methylation ratio matrix (bins × samples), TSV (enabled by `--save-bin-matrix`).

Figure:
- `<out-prefix>.corr.triangle.<pdf|png>`  
  Triangular correlation figure (single or multi-panel).

---

## 3) Installation / requirements

Python packages:
- `numpy`
- `pandas`
- `matplotlib`

Typical install:

```bash
pip install numpy pandas matplotlib
```

---

## 4) Usage

### 4.1 Single panel (one context)

Compute correlation across three samples (e.g., CG):

```bash
python 14_methylation_bin_correlation_triangle.py   --indir CX   --files "S2_CG.cx,S4_CG.cx,S6_CG.cx"   --names "S2,S4,S6"   --bin-size 10000   --min-bin-depth 10   --method pearson   --out-prefix results/CG_10kb
```

This produces:
- `results/CG_10kb.group1.corr.tsv`
- `results/CG_10kb.corr.triangle.pdf`

### 4.2 Multi-panel (CG | CHG | CHH stacked)

```bash
python 14_methylation_bin_correlation_triangle.py   --indir CX   --files "S2_CG.cx,S4_CG.cx,S6_CG.cx|S2_CHG.cx,S4_CHG.cx,S6_CHG.cx|S2_CHH.cx,S4_CHH.cx,S6_CHH.cx"   --group-names "CG|CHG|CHH"   --names "S2,S4,S6"   --bin-size 10000   --min-bin-depth 10   --method pearson   --cmap viridis   --vrange percentile --vpercentile "5,95"   --out-prefix results/CG_CHG_CHH_10kb
```

### 4.3 Save the per-bin methylation ratio matrix

```bash
python 14_methylation_bin_correlation_triangle.py   --indir CX   --files "S2_CG.cx,S4_CG.cx,S6_CG.cx"   --bin-size 20000   --min-bin-depth 20   --save-bin-matrix   --out-prefix results/CG_20kb_depth20
```

### 4.4 Switch triangle and add annotations

```bash
python 14_methylation_bin_correlation_triangle.py   --indir CX   --files "S2_CG.cx,S4_CG.cx,S6_CG.cx"   --triangle upper   --annot   --out-prefix results/CG_upper_annot
```

---

## 5) Arguments (detailed)

### Required

- `--indir`  
  Input directory containing CX files.

- `--files`  
  Files to correlate.

  Syntax:
  - Use `,` to separate files **within** one panel.
  - Use `|` to separate **panels** (multi-panel plot).

  Example:
  - One panel: `A.cx,B.cx,C.cx`
  - Three panels: `A_CG.cx,B_CG.cx|A_CHG.cx,B_CHG.cx|A_CHH.cx,B_CHH.cx`

- `--out-prefix`  
  Output prefix, e.g. `results/CG_10kb`.

### Optional panel labels

- `--group-names`  
  Panel titles separated by `|` (must match panel count).  
  Example: `CG|CHG|CHH`.

- `--names`  
  Comma-separated sample labels for both slanted edges.  
  If not set, the script uses file basenames.

### Binning & filters

- `--bin-size` (default: `10000`)  
  Bin size in bp.

- `--min-bin-depth` (default: `10`)  
  For each (chrom, bin), sum all depths across cytosines; drop bins with summed depth < threshold.
  Increase this to reduce noise in low-coverage bins.

- `--chunksize` (default: `2000000`)  
  Read input in chunks to reduce memory usage.

### Correlation

- `--method {pearson,spearman,kendall}` (default: `pearson`)  
  Correlation method.

### Plotting

- `--fig-format {pdf,png}` (default: `pdf`)  
- `--dpi` (default: `300`)  
  Only used when saving PNG.

- `--cmap` (default: `viridis`)  
  Matplotlib colormap.

- `--triangle {lower,upper}` (default: `lower`)  
  Which triangle to draw.

- `--annot`  
  Annotate each diamond cell with the correlation value.

### Color scaling (vmin / vmax)

The script chooses vmin/vmax from **off-diagonal** correlation values.

- `--vrange {fixed,data,percentile}` (default: `percentile`)
  - `fixed`: use [-1, 1]
  - `data`: use min/max of off-diagonal values
  - `percentile`: use percentiles of off-diagonal values (recommended)

- `--vpercentile` (default: `5,95`)  
  Percentiles used when `--vrange percentile`.

- `--vmin`, `--vmax`  
  Manual limits. **Both must be set** to take effect.

---

## 6) Interpretation notes

- The plot uses **diamond** polygons; only one triangle is drawn.
- The diagonal is not emphasized; correlation self-self is not the focus for comparison.
- Because bins are intersected across all samples, the correlation reflects *shared covered bins*.

---

## 7) Troubleshooting

### “No common (chrom,bin) across samples”
The script will error if the inner-joined bin matrix has 0 rows.

Try:
- Lower `--min-bin-depth` (e.g., 10 → 1 or 5)
- Increase `--bin-size` (e.g., 10kb → 20kb/50kb)
- Verify all files share the same chromosome naming (e.g., `Chr01` vs `chr1`).

### Plot is too “saturated”
- Use percentile scaling:
  - `--vrange percentile --vpercentile 5,95` (default)
- Or manual:
  - `--vmin 0.80 --vmax 1.00` (remember: set both)

### Very large files / slow run
- Increase `--chunksize` (within RAM limits)
- Use larger bins (`--bin-size`) to reduce total bins

---

## 8) Best practices

- Keep `--names` consistent across panels when plotting CG/CHG/CHH together.
- Choose `--bin-size` with genome scale and coverage in mind:
  - High-coverage WGBS: 5–20 kb
  - Low-coverage: 20–100 kb (and reduce `--min-bin-depth`)

