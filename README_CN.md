# OmicsCanvas（中文总览）

OmicsCanvas 是一个按“步骤编号”组织的基因中心多组学分析与绘图流程（ChIP-seq/ATAC-seq/WGBS 甲基化 + RNA 表达）。
共 16 个步骤：**1–4 为准备步骤**，**5–9 为计算步骤**，**10–16 为绘图/分析步骤**。

## 安装

### 方案A：解压即用（推荐）

```bash
pip install -r requirements.txt
python omicscanvas.py --list
python omicscanvas.py 10 -- --help
```

### 方案B：pip 安装（用于发布到 PyPI / GitHub）

```bash
pip install .
omicscanvas --list
omicscanvas run 10 -- --help
```

## 流程总览

| 步骤 | 阶段 | 目的 | 脚本 |
|---:|---|---|---|
| 01 | 准备 | GFF3/GTF -> BED + gene length table | `scripts/01_prepare_gff_to_bed_genes_length.py` |
| 02 | 准备 | Split CX report by context (CG/CHG/CHH) | `scripts/02_prepare_cx_context_split.py` |
| 03 | 准备 | Merge CX replicates | `scripts/03_prepare_cx_replicate_merge.py` |
| 04 | 准备 | Extract per-gene methylation table from CX | `scripts/04_prepare_extract_gene_methylation.py` |
| 05 | 计算 | BAM -> gene matrix (TSS/gene/TES) | `scripts/05_compute_bam_to_gene_matrices.py` |
| 06 | 计算 | Build CX_gene matrices/profiles from CX/genes | `scripts/06_compute_cx_gene_matrix.py` |
| 07 | 计算 | BAM -> FPKM | `scripts/07_compute_bam_to_fpkm.py` |
| 08 | 计算 | Differential expression (NB/other) helper | `scripts/08_compute_nb_de.py` |
| 09 | 计算 | Merge expression + pairwise DE pipeline | `scripts/09_compute_merge_expr_pairwise_de.py` |
| 10 | 绘图 | Plot whole profiles in 2D/3D | `scripts/10_plot_whole_profile_2d3d.py` |
| 11 | 绘图 | Track vs expression heatmap | `scripts/11_plot_track_vs_expr_heatmap.py` |
| 12 | 绘图 | Histone clustering pipeline from matrices | `scripts/12_plot_histone_cluster_pipeline.py` |
| 13 | 绘图 | Methylation profile 2D/3D | `scripts/13_plot_methylation_profile_2d3d.py` |
| 14 | 绘图 | Methylation vs expression heatmap | `scripts/14_plot_meth_vs_expr_heatmap.py` |
| 15 | 绘图 | Gene tracks (2D/3D) | `scripts/15_plot_gene_tracks_2d3d.py` |
| 16 | 绘图 | Gene circle plot | `scripts/16_plot_gene_circle_plot.py` |

> 详细英文说明已合并在 `README_EN.md`（内容更完整）；每个脚本也支持 `--help` 查看参数。
