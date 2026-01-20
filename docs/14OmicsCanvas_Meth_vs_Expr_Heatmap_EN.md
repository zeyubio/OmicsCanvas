# OmicsCanvas — Methylation vs Expression Heatmap

This tool integrates:
- **Methylation segment tables** (upstream / body / downstream)
- **Expression table (FPKM/TPM)**

and produces:
1) A **heatmap** showing methylation ratio along the merged profile after ranking genes by expression
2) A saved PDF figure for publication / review

---

## 1) Inputs

### 1.1 Methylation files (3 segments)
The script searches for three files:

```
<meth_dir>/<sample>_<cx><suffix-start>
<meth_dir>/<sample>_<cx><suffix-gene>
<meth_dir>/<sample>_<cx><suffix-end>
```

Example (your current naming):
- `CX_gene/SRR9321764_CHH_upstream_bins50.tsv`
- `CX_gene/SRR9321764_CHH_body_bins100.tsv`
- `CX_gene/SRR9321764_CHH_downstream_bins50.tsv`

Use:
- `--meth-dir CX_gene`
- `--suffix-start _upstream_bins50.tsv`
- `--suffix-gene  _body_bins100.tsv`
- `--suffix-end   _downstream_bins50.tsv`

#### Required columns
Each segment file must contain columns:

| column | meaning |
|---|---|
| `name` | gene ID |
| `po`   | bin index (can be any continuous integer range; the script will normalize by subtracting min) |
| `me`   | methylated counts |
| `al`   | total counts |

If the same `(name, po)` appears multiple times, the script sums `me` and `al` before computing ratios.

#### Segment length rules
- By default, the script requires the three segments to have the **same number of bins**.
- If your segments have different lengths (e.g., 50/100/50 is fine if each segment file itself has a consistent range), enable:

`--allow-unequal-segments`

The boundary positions (TSS/TES) will follow the inferred segment lengths.

### 1.2 Expression table
A TSV file where:
- **Row index (first column)** is gene ID
- Columns are replicates or samples

Use:
- `--fpkm FPKM.txt`
- `--fpkm-cols 0,1,2` (0-based column indices to average)

---

## 2) What the script does

1) Load upstream/body/downstream methylation tables
2) Normalize each segment `po` to start from 0 (subtract segment `po_min`)
3) Offset segments and concatenate them into one global coordinate
4) Compute methylation ratio per `(gene, global_po)`:

`ratio = me / al`

5) Load expression, compute mean across `--fpkm-cols`
6) Rank genes by expression:
   - genes with expression <= `--zero-threshold` go into `--none-bins` bins
   - genes with expression >  `--zero-threshold` go into `--exp-bins` bins
7) For each expression bin, aggregate methylation ratio along global bins
8) Plot a heatmap

---

## 3) Outputs

The figure is written to:

`<out_prefix>_<sample>_<cx>_meth_vs_expr_heatmap.pdf`

Saved in the current working directory.

---

## 4) Example command (your current naming)

```bash
python omicscanvas_meth_vs_expr_heatmap_EN.py \
  --sample SRR9321764 \
  --cx CHH \
  --meth-dir CX_gene \
  --suffix-start _upstream_bins50.tsv \
  --suffix-gene  _body_bins100.tsv \
  --suffix-end   _downstream_bins50.tsv \
  --fpkm FPKM.txt \
  --fpkm-cols 0,1,2 \
  --none-bins 10 \
  --exp-bins 90 \
  --distance 2000 \
  --scale-mode ratio \
  --cmap RdBu_r \
  --vmax 0.25 \
  --out-prefix SRR9321764
```

---

## 5) Parameter reference (common ones)

- `--sample` : sample ID for filename prefix.
- `--cx` : context (CG/CHG/CHH).
- `--meth-dir` : directory for methylation segment tables.
- `--suffix-start/--suffix-gene/--suffix-end` : file suffix for each segment.
- `--fpkm` : expression table.
- `--fpkm-cols` : 0-based indices to average.
- `--none-bins` : number of bins for genes with expression <= threshold.
- `--exp-bins` : number of bins for genes with expression > threshold.
- `--zero-threshold` : threshold (default 0.0).
- `--distance` : label value shown at profile edges (e.g. -2000bp / +2000bp).
- `--scale-mode` :
  - `ratio` (default): use raw ratio, typically [0, 1]
  - `diverging`: symmetric scale around 0 (useful after transforming)
  - `quantile`: auto-range from quantiles
- `--vmin/--vmax` : manually set color scale bounds.
- `--range` : convenience bound (if given, sets both vmin/vmax depending on mode).
- `--cmap` : colormap.
- `--out-prefix` : output prefix.

---

## 6) Troubleshooting

- **Heatmap is empty/blank**: most often caused by no overlap between methylation gene IDs and expression gene IDs.
- **Division by zero**: bins with `al=0` become NaN; you can fill with `--fillna 0`.
- **Segment bins error**: if segment lengths differ, enable `--allow-unequal-segments`.

