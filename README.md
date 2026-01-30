# OmicsCanvas: A Multi-Omics Integration & Visualization Toolkit

[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](CONTRIBUTING.md)

> **A unified framework for integrating and visualizing epigenetic, epitranscriptomic, and transcriptomic data in a gene-centric coordinate system.**

## 📖 Introduction

**OmicsCanvas** is a one-stop toolkit designed to bridge the gap between gene expression and diverse regulatory layers (e.g., ChIP-seq, ATAC-seq, WGBS, m6A-seq). Unlike traditional peak-centric tools, OmicsCanvas employs a **unified gene-centric coordinate framework**, decomposing genes into **promoter, body, and flanking neighborhoods**.

Key capabilities include:
* **Unified Binning:** Aligns heterogeneous signals (BAM/BigWig) and base-resolution methylation data (Bismark CX) to a shared structure.
* **Expression Integration:** Directly links epigenomic signals with gene expression levels (TPM/FPKM).
* **Publication-Ready Visualization:** Exports high-quality **2D, Pseudo-3D, and Circular** vector graphics.
* **Unsupervised Learning:** Identifies cooperative or antagonistic regulatory patterns via K-means clustering.

---

## 📦 Installation

OmicsCanvas follows a minimal-dependency philosophy.

### Option A: Conda / Mamba (Recommended)
```bash
# Create environment
mamba create -n omicscanvas python=3.9 -y
mamba activate omicscanvas

# Install core scientific stack & bioinformatics tools
mamba install -c conda-forge numpy pandas matplotlib scipy scikit-learn -y
mamba install -c bioconda pysam -y
```

### Option B: Pip

```bash
pip install numpy pandas matplotlib scipy scikit-learn pysam
```


## 🚀 Usage & Workflow
The OmicsCanvas workflow consists of three main stages: Preparation, Matrix Calculation, and Visualization.


### Step 1: Preparation
Convert your annotation (GFF3) into standard BED format and calculate gene lengths.
```bash

python scripts/omicscanvas_gff_to_bed_genes_length.py 
  -i annotation.gff3 \
  --bed-feature gene \
  --length-feature CDS \
  --bed-zero-based \
  -o gene.bed \
  -l gene_cds_length.tsv \
  --merge-overlaps
```

### Step 2: Matrix Generation
A. For BAM Files (ChIP / ATAC / RNA)
Convert sorted BAM files into gene-centric coverage matrices (Promoter, Body, Terminator).

```bash

python scripts/omicscanvas_bam_to_gene_matrices.py \
  -b sample.sorted.bam \
  -g gene.bed \
  --outdir matrices \
  --distance 2000 \
  --tss-bins 100 --gene-body-bins 100 --tes-bins 100 \
  -o Sample_H3K4me3
```

> **Output:** This process generates three core matrix files: `_tss_matrix.tsv`, `_gene_profile_matrix.tsv`, and `_tes_matrix.tsv`.


### B. For Methylation (WGBS)
Process Bismark CX reports into gene-centric methylation matrices.

```bash

# 1. Split CX report by context
python scripts/cx_context_split.py -i sample.CX_report.txt -p sample -d meth_data

# 2. Generate gene matrix (e.g., CHH context)
python scripts/cx_gene_matrix.py -s sample -c CHH -b gene.bed --cx-dir meth_data -o CX_matrices
```


> **Output:** This process generates three core matrix files: `_tss_matrix.tsv`, `_gene_profile_matrix.tsv`, and `_tes_matrix.tsv`.

### 🎨 Visualization Gallery
1. Global Multi-Omics Profile (Pseudo-3D)
Visualize the genome-wide distribution of histone modifications or accessibility. The 3D mode allows stacking multiple tracks for intuitive comparison.

```bash

python scripts/omicscanvas_plot_whole_profile_2d3d.py \
  --mode 3d \
  --matrix-dir matrices \
  --gene-type gene \
  --group "G_H3K4me3,B_H3K4me3|G_H3K27me3,B_H3K27me3" \
  --ylabels "H3K4me3|H3K27me3" \
  --out results/global_profile_3d.png
 ```

<div align="center">
  <img src="./images/fig1_global_3D_1.png" width="800px">
  <p><b>Figure 1:</b> Global Pseudo-3D profile showing multi-omics signal distribution across the genome.</p>
</div>


### 2. Signal vs. Expression Heatmap
Sort genes by expression level (High to Low) and visualize the corresponding epigenetic signal density.


```bash

python scripts/omicscanvas_histone_vs_expr_heatmap.py \
  --matrix-dir matrices \
  --tracks "H3K4me3" \
  --expr expression_FPKM.txt \
  --exp-bins 90 \
  --cmap RdBu_r \
  --out-prefix results/H3K4me3_vs_Expr
 ```

<table style="width: 100%; text-align: center; border-collapse: collapse; border: none;">
  <tr>
    <td style="border: none; width: 33%;">
      <img src="./images/fig2_heatmap_H3K4me3_1.png" width="100%">
      <br>
      <p><i>Sample 1: H3K4me3 TSS Signal</i></p>
    </td>
    <td style="border: none; width: 33%;">
      <img src="./images/fig2_heatmap_H3K4me3_2.png" width="100%">
      <br>
      <p><i>Sample 2: H3K4me3 GeneBody Signal</i></p>
    </td>
    <td style="border: none; width: 33%;">
      <img src="./images/fig2_heatmap_H3K4me3_3.png" width="100%">
      <br>
      <p><i>Sample 3: H3K4me3 TES Signal</i></p>
    </td>
  </tr>
</table>


### 3. Clustering Analysis
Use K-means clustering to identify distinct chromatin states or regulatory patterns across samples.

```bash

python scripts/omicscanvas_histone_cluster_pipeline.py \
  --matrix-dir matrices \
  --in-group "Sample_H3K4me3;Sample_H3K27me3" \
  --k 4 \
  --out-prefix results/Clustering_Analysis
```


<table style="width: 100%; text-align: center; border-collapse: collapse; border: none;">
  <tr>
    <td style="border: none; width: 33%;">
      <img src="./images/fig3_cluster_heatmap_TSS.png" width="100%">
      <br>
      <p><i>Sample 1: H3K4me1 H3K4me3 TSS cluster heatmap Signal</i></p>
    </td>
    <td style="border: none; width: 33%;">
      <img src="./images/fig3_cluster_heatmap_genebody.png" width="100%">
      <br>
      <p><i>Sample 2: H3K4me1 H3K4me3 Genebody cluster heatmap Signal</i></p>
    </td>
    <td style="border: none; width: 33%;">
      <img src="./images/fig3_cluster_heatmap_TES.png" width="100%">
      <br>
      <p><i>Sample 3: H3K4me1 H3K4me3 TES cluster heatmap Signal</i></p>
    </td>
  </tr>
</table>


<table style="width: 100%; text-align: center; border-collapse: collapse; border: none;">
  <tr>
    <td style="border: none; width: 33%;">
      <img src="./images/fig3_cluster_lineplot_TSS.png" width="100%">
      <br>
      <p><i>Sample 1: H3K4me1 H3K4me3 TSS cluster heatmap Signal</i></p>
    </td>
    <td style="border: none; width: 33%;">
      <img src="./images/fig3_cluster_lineplot_genebody.png" width="100%">
      <br>
      <p><i>Sample 2: H3K4me1 H3K4me3 Genebody cluster heatmap Signal</i></p>
    </td>
    <td style="border: none; width: 33%;">
      <img src="./images/fig3_cluster_lineplot_TES.png" width="100%">
      <br>
      <p><i>Sample 3: H3K4me1 H3K4me3 TES cluster heatmap Signal</i></p>
    </td>
  </tr>
</table>


<div align="center">
  <img src="./images/fig3_cluster_boxplot.png" width="600px">
  <p><b>Figure: Gene Expression Distribution Across Clusters</b><br>
  <i>Boxplots showing the TPM/FPKM expression levels for each identified K-means cluster, highlighting the correlation between epigenetic signals and transcriptional activity.</i></p>
</div>


### 4. Single Gene Visualization (Pseudo-3D)
Zoom in on specific candidate genes. Stack ChIP-seq and RNA-seq tracks in a 3D layout to show co-occupancy.

```bash
python scripts/omicscanvas_gene_tracks_2d3d.py \
  --mode 3d \
  --gff3 annotation.gff3 \
  --gene Potri.006G061800 \
  --bam-dir bam_files \
  --bam-spec "H3K4me3.bam;H3K27me3.bam" \
  --out results/candidate_gene_3d.png
```

<div align="center">

  <h3>📍 Single-Gene Multi-Omics 2D Track</h3>
  <img src="./images/fig4_gene_track_2D .png" width="750px" alt="2D Gene Track">
  <p><i>A high-resolution 2D visualization showing the distribution of epigenetic signals across the gene body and regulatory elements.</i></p>

  <br>

  <h3>🧊 Single-Gene Multi-Omics 3D Track</h3>
  <img src="./images/fig4_gene_track_3D .png" width="750px" alt="3D Gene Track">
  <p><i>Enhanced Pseudo-3D perspective allowing for intuitive comparison of stacked signal intensities and co-occupancy patterns.</i></p>

</div>


### 5. Circular Gene Plot
Map genomic windows onto angular coordinates. Ideal for visualizing complex multi-layer regulation (e.g., Histone + Methylation) in a compact format.

```bash
python scripts/omicscanvas_gene_circle_plot.py \
  --gff3 annotation.gff3 \
  --gene Potri.006G061800 \
  --bam-spec "H3K4me3.bam;ATAC.bam" \
  --meth-spec "sample_CHH" \
  --circle-meth-layout combined \
  --out results/candidate_gene_circle.png
```

<div align="center">
  <h2>⭕ Interactive Gene Circular Visualization</h2>
  <img src="./images/fig5_circle_gene_track_3D .png" width="600px" alt="Circular Gene Track">
  <br>
  <p align="center" style="width: 80%;">
    <b>Figure: Integrative Circular Track of Histone Modifications and Transcriptomics</b><br>
    <i>This circular coordinate framework maps multi-layer regulatory data (e.g., ChIP-seq signals and RNA-seq expression) onto an angular axis. It provides a compact yet comprehensive view of the epigenetic landscape and transcriptional activity for a specific candidate gene.</i>
  </p>
</div>



### 📚 Citation
If you use OmicsCanvas in your research, please cite:

OmicsCanvas: A multi-omics platform for integration and visualization of epigenetic and epitranscriptomic regulation. Zeyu Zhang, Zuoling Ma, Tian Hua, et al.

### ✉️ Contact
Issues: Please post feature requests or bugs to the Issues page.

Email: lfgu@fafu.edu.cn
