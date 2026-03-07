# 16_plot_meth_vs_expr_heatmap.py — Methylation vs Expression Heatmap

This script ranks genes by expression (FPKM/TPM/count-like), splits genes into expression bins, then calculates **methylation ratio** for each bin along a **merged gene meta-profile** (upstream / gene body / downstream) and plots a heatmap.

- Methylation ratio per bin position: **ratio = sum(me) / sum(al)** (computed across genes within each expression bin)
- X-axis: upstream → **TSS** → **TES** → downstream (segment boundaries are derived from your methylation segment tables)
- Y-axis: expression bins (low → high)

---

## 0) Quick sanity check: the README you provided is for a different script

Your current README `16OmicsCanvas_StepX_Track_vs_Expr_Heatmap_EN.updated.md` describes **Step2 track matrices** and explicitly names `12_plot_histone_vs_expr_heatmap.py` plus `--matrix-dir/--tracks/--gene-types` etc. This does **not** match this methylation script.  
See: it states the script file is `12_plot_histone_vs_expr_heatmap.py` and explains Step2 matrices + `--tracks` syntax. fileciteturn19file0L1-L13 fileciteturn19file0L57-L82

This new README below matches the actual CLI options and file conventions implemented in `16_plot_meth_vs_expr_heatmap.py`. fileciteturn18file0L445-L610

---

## 1) Inputs

### 1.1 Methylation segment tables (required)

The script reads **three tables** and merges them into one global position axis:

- `<meth-dir>/<sample>_<cx><suffix-start>`  (upstream/TSS segment)
- `<meth-dir>/<sample>_<cx><suffix-gene>`   (gene-body segment)
- `<meth-dir>/<sample>_<cx><suffix-end>`    (downstream/TES segment)

Default suffixes are:
- `--suffix-start _gene_start_300.txt`
- `--suffix-gene  _gene_gene_300.txt`
- `--suffix-end   _gene_end_300.txt` fileciteturn18file0L466-L485

So, by default, the expected filenames look like:
- `gene/<sample>_<cx>_gene_start_300.txt`
- `gene/<sample>_<cx>_gene_gene_300.txt`
- `gene/<sample>_<cx>_gene_end_300.txt` fileciteturn18file0L164-L171

**Required columns (TSV):**
- `bin` : position/bin index inside the segment (can start at 0/1/any integer; the script normalizes within segment)
- `me`  : methylated counts
- `al`  : total coverage/depth
- `name`: gene ID (can be either a column OR stored in the first column as the dataframe index)

Important details:
- The script reads with `index_col=0` first; if there is no `name` column, it uses the index as gene ID. fileciteturn18file0L100-L107
- The script **requires a column named `bin`** (it then converts it to `po`). fileciteturn18file0L108-L120

> Note: the docstring mentions `po`, but the current implementation expects the input column to be called `bin`.

### 1.2 Expression table (required)

- TSV file with **gene IDs as row index** (first column)
- Choose which columns to average with `--fpkm-cols` (0-based indices, **excluding** the index column). fileciteturn18file0L495-L507

Example:
```
GeneID	rep1	rep2	rep3
GeneA	 1.2	 1.1	 1.3
GeneB	 0.0	 0.0	 0.0
...
```

The script intersects genes between methylation tables and expression table and fails if there is no overlap. fileciteturn18file0L580-L586

---

## 2) What the script does

### 2.1 Merge three methylation segments into one axis

For each segment (start/gene/end), the script:
1) finds `po_min/po_max` from `bin`
2) normalizes positions to `0..seg_bins-1`
3) offsets gene/end segments cumulatively so all positions are on one global axis
4) sums duplicated (gene,position) entries by `(name, po)`.

By default it requires equal bin counts across the three segments, otherwise it errors; you can override this with `--allow-unequal-segments`. fileciteturn18file0L177-L184 fileciteturn18file0L486-L493

### 2.2 Bin genes by expression

Genes are sorted from low → high expression.
- `expr <= --zero-threshold` → split into `--none-bins` bins
- `expr >  --zero-threshold` → split into `--exp-bins` bins fileciteturn18file0L510-L513 fileciteturn18file0L234-L257

Each bin produces **one heatmap row**.

### 2.3 Build the heatmap matrix

For each expression bin and each global position:
- `ratio = sum(me) / sum(al)` across all genes in that bin
- If `al==0`, ratio is `NaN`; optionally fill NaNs with `--fillna`. fileciteturn18file0L279-L290 fileciteturn18file0L548-L553

---

## 3) Outputs

The output filename is always a **PDF**:

`<out-prefix>_<sample>_<cx>_meth_vs_expr_heatmap.pdf` fileciteturn18file0L608-L624

- If `--out-prefix` is not provided, it defaults to `--sample`. fileciteturn18file0L557-L610

---

## 4) Installation / requirements

Python packages:
- numpy
- pandas
- seaborn
- matplotlib fileciteturn18file0L52-L56

Install:
```bash
pip install numpy pandas seaborn matplotlib
```

---

## 5) Usage

### 5.1 Minimal example (default suffixes)

```bash
python 16_plot_meth_vs_expr_heatmap.py   --sample SRR9321764   --cx CHH   --meth-dir gene   --fpkm FPKM.txt   --fpkm-cols 0,1,2   --out-prefix SRR9321764
```

### 5.2 Using your “new naming” suffixes (recommended)

If your methylation tables are named like:
- `CX_gene/SRR9321764_CHH_upstream_bins50.tsv`
- `CX_gene/SRR9321764_CHH_body_bins100.tsv`
- `CX_gene/SRR9321764_CHH_downstream_bins50.tsv`

Run:

```bash
python 16_plot_meth_vs_expr_heatmap.py   --sample SRR9321764   --cx CHH   --meth-dir CX_gene   --suffix-start _upstream_bins50.tsv   --suffix-gene  _body_bins100.tsv   --suffix-end   _downstream_bins50.tsv   --fpkm FPKM.txt   --fpkm-cols 0,1,2   --none-bins 10   --exp-bins 90   --zero-threshold 0   --scale-mode ratio --range 1   --cmap RdBu_r   --out-prefix CHH_bins50_100_50
```

Suffix handling and naming notes are in the script help and description. fileciteturn18file0L466-L485 fileciteturn18file0L445-L448

---

## 6) Arguments (complete)

### Required

- `--sample`  
  Sample ID used to construct filenames: `<meth-dir>/<sample>_<cx><suffix*>`. fileciteturn18file0L450-L455
- `--cx`  
  Methylation context: `CG`, `CHG`, or `CHH`. fileciteturn18file0L457-L460
- `--fpkm`  
  Expression table path (TSV). fileciteturn18file0L495-L501

### File layout

- `--meth-dir` (default: `gene`)  
  Directory containing methylation segment tables. fileciteturn18file0L461-L464
- `--suffix-start` / `--suffix-gene` / `--suffix-end`  
  Suffixes for start/gene/end tables. fileciteturn18file0L466-L485
- `--allow-unequal-segments`  
  Allow different bin counts across start/gene/end. fileciteturn18file0L486-L493

### Expression columns & binning

- `--fpkm-cols` (default: `0,1,2`)  
  0-based columns (excluding index column) used to compute mean expression. fileciteturn18file0L502-L507
- `--none-bins` (default: 10)  
  Number of bins for genes with expression `<= --zero-threshold`. fileciteturn18file0L510-L513
- `--exp-bins` (default: 90)  
  Number of bins for genes with expression `> --zero-threshold`. fileciteturn18file0L510-L513
- `--zero-threshold` (default: 0)  
  Threshold for none/expressed split. fileciteturn18file0L510-L513

### Plot axis labels

- `--distance` (default: 2000)  
  Only affects x-axis labels (`-2000bp`, `TSS`, `TES`, `+2000bp`); does not affect bin merging. fileciteturn18file0L514-L519 fileciteturn18file0L414-L418

### Color scaling

- `--scale-mode {ratio,diverging,quantile}` (default: `ratio`) fileciteturn18file0L522-L527
- `--range FLOAT`  
  Shortcut range:
  - `ratio/quantile`: `[0, range]`
  - `diverging`: `[-range, +range]` fileciteturn18file0L529-L534
- `--vmin`, `--vmax`  
  Manual limits (use both together). fileciteturn18file0L536-L538
- `--quantiles "q_low,q_high"` (default: `0.01,0.99`)  
  Used only in `scale-mode=quantile`. fileciteturn18file0L537-L539
- `--no-clamp1`  
  In `ratio` mode, do not clamp `vmax` to `<= 1`. fileciteturn18file0L538-L540
- `--cmap`  
  Colormap name. If not set, defaults to `viridis` (ratio/quantile) or `RdBu_r` (diverging). fileciteturn18file0L297-L301 fileciteturn18file0L541-L546

### Missing values & cosmetics

- `--fillna FLOAT`  
  Fill NaNs with a constant. fileciteturn18file0L548-L553
- `--no-box`  
  Disable outer border box. fileciteturn18file0L554-L555
- `--ytick-step INT`  
  Control y tick spacing (default auto). fileciteturn18file0L554-L556
- `--out-prefix`  
  Output prefix (default: `--sample`). fileciteturn18file0L557-L610

---

## 7) Troubleshooting

### “missing columns: ['bin', ...]”
Your methylation tables must contain `bin, me, al` (and `name` either as a column or as row index). fileciteturn18file0L108-L120  
If your input uses `po` instead of `bin`, rename the column to `bin`.

### “Segment bins are not equal”
By default the script requires equal bin counts across start/gene/end segments; use `--allow-unequal-segments` to proceed. fileciteturn18file0L177-L184

### “No overlap genes between methylation and expression table”
Gene IDs must match exactly; this script does not implement `strip_dot`. The overlap check is strict. fileciteturn18file0L580-L586

### Heatmap has blank cells
Those positions have `al==0` (or missing bins). Consider:
- lowering `--zero-threshold` or adjusting bins
- setting `--fillna 0`
- increasing coverage / filtering input tables

---

## 8) Notes for reproducibility

- Keep `--none-bins + --exp-bins` constant across runs if you want comparable y-axis resolution.
- In `scale-mode=ratio`, the script clamps `vmax <= 1` by default; use `--no-clamp1` if you expect ratios > 1 due to upstream preprocessing.
