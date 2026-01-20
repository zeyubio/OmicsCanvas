# OmicsCanvas — Methylation whole-profile plotting (2D / 3D)

This tool plots **methylation “whole-profile” curves** (typically CG/CHG/CHH) along a gene meta-profile, in either:
- **2D**: all contexts stacked vertically in one figure, each context showing multiple samples as lines
- **3D**: layered (pseudo-3D) stacking of contexts to visually separate layers

The script reads per-sample, per-context profile files and produces PDF figures.

---

## 1) Input files

### 1.1 File naming rule
For each sample and context, the script looks for:

```
<meth_dir>/<sample>_<context><whole_suffix>
```

Examples:
- `CX_gene/SRR9321764_CHH_profile.tsv`  (use `--meth-dir CX_gene --whole-suffix _profile.tsv`)
- `gene/LSH10_CG_whole_line.txt`        (use `--whole-suffix _whole_line.txt`)

### 1.2 File format
Each profile file should contain **one numeric column** (no header required). The script loads it as a single column named `methylation`.

The expected profile length is:

```
all_length = bins_st + bins_body + bins_en
```

If your profile length is different, set `--bins-st/--bins-body/--bins-en` to match how you generated the profile.

---

## 2) Outputs

Depending on `--mode`, the script writes:
- `<out_prefix>.2D.pdf`
- `<out_prefix>.3D.pdf`

Files are saved to the **current working directory**.

---

## 3) Common usage examples

### 3.1 Plot CHH profile for a single sample (2D)

```bash
python omicscanvas_methylation_profile_2d3d_EN.py \
  --meth-dir CX_gene \
  --contexts CHH \
  --samples SRR9321764 \
  --labels SRR9321764 \
  --whole-suffix _profile.tsv \
  --mode 2d \
  --out-prefix SRR9321764_CHH_profile
```

### 3.2 Plot CG/CHG/CHH for one sample (3D)

```bash
python omicscanvas_methylation_profile_2d3d_EN.py \
  --meth-dir CX_gene \
  --contexts CG,CHG,CHH \
  --samples SRR9321764 \
  --labels SRR9321764 \
  --whole-suffix _profile.tsv \
  --mode 3d \
  --out-prefix SRR9321764_CX_profile_3D
```

### 3.3 Plot multiple samples per context

```bash
python omicscanvas_methylation_profile_2d3d_EN.py \
  --meth-dir CX_gene \
  --contexts CHH \
  --samples SRR9321764,SRR9321765 \
  --labels WT,mut \
  --whole-suffix _profile.tsv \
  --mode both \
  --out-prefix CHH_WT_vs_mut
```

### 3.4 Provide manual line colors

```bash
python omicscanvas_methylation_profile_2d3d_EN.py \
  --meth-dir CX_gene \
  --contexts CHH \
  --samples SRR9321764,SRR9321765 \
  --labels WT,mut \
  --whole-suffix _profile.tsv \
  --line-colors '#1f77b4,#d62728' \
  --out-prefix CHH_color_demo
```

---

## 4) Parameter reference (public/common parameters)

### Input / sample definition
- `--meth-dir`  
  Directory that stores profile files.

- `--contexts`  
  Comma-separated contexts (e.g. `CG,CHG,CHH`). Each context becomes one layer/row.

- `--samples`  
  Comma-separated sample IDs used in the filename prefix `<sample>_<context>...`.

- `--labels`  
  Comma-separated legend labels for samples. Must match `--samples` length.

- `--whole-suffix`  
  Suffix appended after `<sample>_<context>`.  
  Examples: `_profile.tsv`, `_whole_line.txt`.

### Gene meta-profile layout (used for ticks/vertical lines)
- `--distance`  
  Only used for x-axis labeling (e.g. `-2000bp`/`+2000bp`). Not used for scaling.

- `--bins-st`  
  Number of bins in upstream (before TSS).

- `--bins-body`  
  Number of bins in gene body (TSS → TES).

- `--bins-en`  
  Number of bins in downstream (after TES).

### Figure & style
- `--fig-x`, `--fig-y`  
  Matplotlib figure size (inches).

- `--cmap`  
  Colormap name used to auto-assign line colors when `--line-colors` is not given.

- `--line-colors`  
  Comma-separated colors for samples, e.g. `'#1f77b4,#d62728'`. Colors cycle if fewer than samples.

- `--line-lw`  
  Line width.

- `--ylim`  
  Y-axis range(s). Either a single `min,max` for all contexts, or semicolon-separated per context.  
  Example: `0,0.3` or `0,0.3;0,0.2;0,0.1`.

- `--mode`  
  `2d`, `3d`, or `both`.

### 3D-only layout controls
- `--edge`  
  Outer margin (0–0.5 recommended). Larger = more whitespace.

- `--x-offset`, `--y-offset`  
  Layer offsets in percent (0–100). Controls pseudo-3D stacking.

### Output
- `--out-prefix`  
  Output file prefix.

- `--standard-distance-labels`  
  If set, TSS/TES labels are shown as `+distance ... -distance` style.

---

## 5) Troubleshooting

- **No output file appears**: outputs are written to the *current directory*; run `ls -lh *.pdf` after execution.
- **File not found**: confirm the naming rule `<sample>_<context><whole_suffix>` and set `--meth-dir` and `--whole-suffix`.
- **Wrong tick positions**: `--bins-st/--bins-body/--bins-en` must match the profile generation step.

