# OmicsCanvas — Gene Circle Plot (CN/EN Manual)
**Script:** `16_plot_gene_circle_plot.py`  
**CLI (argparse prog):** `omicscanvas_gene_circle_plot.py`  

> This document is a bilingual (Chinese/English) user manual for OmicsCanvas Script 16.  

---

### 1. Overview
`omicscanvas_gene_circle_plot.py` generates a **single-gene circular (ring) plot** within a fixed window (gene ± distance).

**What it can render:**
1) **Gene structure ring (required):** UTR / CDS / gene-body backbone + strand arrow  
2) **BAM coverage rings (optional):** spike-style coverage computed from BAM files  
3) **Methylation rings (optional):** site-level CG/CHG/CHH ratios from TSV; `combined` overlay or `separate` rings  

**Typical use cases:**
- Publication-grade local signals around one gene (ChIP/ATAC/RNA, etc.)
- Overlay DNA methylation (CG/CHG/CHH) in the same figure
- Batch plotting with consistent borders and comparable scales

---

### 2. Requirements
**Python:** 3.8+ recommended  
**Core packages:** `numpy`, `pandas`, `matplotlib`  
**For BAM rings:** `pysam` is required.

Conda install:
```bash
conda install -c conda-forge numpy pandas matplotlib
conda install -c bioconda pysam
```

---

### 3. Input preparation

#### 3.1 GFF3 (required)
**Critical:** `--gene` must be found in GFF3 as the transcript/mRNA ID.

Recommended:
- `mRNA` feature with `ID=...` (default transcript key)
- `exon`/`CDS`/UTR features with `Parent=...` (default linkage)

If your GFF3 uses different keys, set:
- `--gene-keywords` (default: `ID`)
- `--exon-keywords` (default: `Parent`)

---

#### 3.2 BAM (optional)
BAM must be coordinate-sorted and indexed (`.bai`).

You provide:
- `--bam-dir`
- `--bam-spec`
- `--name-spec`

**Normalization:** coverage is normalized to a comparable scale; `--reads-length` and `--seq-type` strongly affect ring heights.

---

#### 3.3 Methylation TSV (optional)
Provide:
- `--meth-dir`
- `--meth-spec`

Default file naming per prefix:
- `<prefix>__CG.tsv`, `<prefix>__CHG.tsv`, `<prefix>__CHH.tsv`

Recommended TSV:
- header with `po`, `me`, `al` (case-insensitive)
- `po` is a 0-based position within the window: `0 <= po < window_len`
- ratio = `me/al` (filter `al==0` upstream)

If TSV includes `name`, rows are filtered by `name == --gene`. Otherwise, all rows are treated as belonging to the target gene.

---

### 4. Sample spec syntax: `;` and `,`
- `;` separates groups
- `,` separates items within a group

The same syntax applies to `--bam-spec`, `--name-spec`, and `--meth-spec`.

**Hard rules:**
- If `--bam-spec` is set, `--name-spec` must be set
- `--name-spec` must match `--bam-spec` exactly (same groups and item counts per group)

---

### 5. Minimal required parameters
Required:
- `--gff3`, `--gene`, `--distance`, `--out`

For BAM rings:
- `--bam-dir`, `--bam-spec`, `--name-spec`, `--reads-length`, `--seq-type`

For methylation:
- `--meth-dir`, `--meth-spec`, `--circle-meth-layout`

Recommended:
- `--bam-group-ylims`, `--circle-lim`, `--circle-frame`

---

### 6. Quick start examples

#### 6.1 Gene structure only
```bash
python 16_plot_gene_circle_plot_ipynbmethod_groupylim_methfix.py \
  --gff3 genome/Ptrichocarpa_210_v3.0.gene.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --out Potri.001G055900.5.v3.0_circle_gene_only.pdf
```

#### 6.2 BAM rings (no methylation)
```bash
python 16_plot_gene_circle_plot_ipynbmethod_groupylim_methfix.py \
  --gff3 genome/Ptrichocarpa_210_v3.0.gene.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "SRR8742373.sorted.bam,SRR8742374.sorted.bam;SRR8742375.sorted.bam,SRR8742376.sorted.bam;SRR8742314.sorted.bam,SRR8742315.sorted.bam" \
  --name-spec "H3K27me3,H3K36me3;H3K56ac,H3K4me3;RNA_1,RNA_2" \
  --reads-length 150 \
  --seq-type paired \
  --out Potri.001G055900.5.v3.0_circle_no_meth.pdf
```

#### 6.3 BAM + methylation (connect to Script 04 outputs)
Assuming Script 04 produced under `single_gene/`:
- `one__Potri.001G055900.5.v3.0__CG.tsv` (and CHG/CHH)

Use:
- `--meth-dir single_gene`
- `--meth-spec one__Potri.001G055900.5.v3.0`

**(A) combined**
```bash
python 16_plot_gene_circle_plot_ipynbmethod_groupylim_methfix.py \
  --gff3 genome/Ptrichocarpa_210_v3.0.gene.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "SRR8742373.sorted.bam,SRR8742374.sorted.bam;SRR8742375.sorted.bam,SRR8742376.sorted.bam;SRR8742314.sorted.bam,SRR8742315.sorted.bam" \
  --name-spec "H3K27me3,H3K36me3;H3K56ac,H3K4me3;RNA_1,RNA_2" \
  --reads-length 150 \
  --seq-type paired \
  --meth-dir single_gene \
  --meth-spec one__Potri.001G055900.5.v3.0 \
  --circle-meth-layout combined \
  --circle-frame \
  --out Potri.001G055900.5.v3.0_circle_with_meth_combined.pdf
```

**(B) separate**
```bash
python 16_plot_gene_circle_plot_ipynbmethod_groupylim_methfix.py \
  --gff3 genome/Ptrichocarpa_210_v3.0.gene.gff3 \
  --gene Potri.001G055900.5.v3.0 \
  --distance 2000 \
  --bam-dir bam \
  --bam-spec "SRR8742373.sorted.bam,SRR8742374.sorted.bam;SRR8742375.sorted.bam,SRR8742376.sorted.bam;SRR8742314.sorted.bam,SRR8742315.sorted.bam" \
  --name-spec "H3K27me3,H3K36me3;H3K56ac,H3K4me3;RNA_1,RNA_2" \
  --reads-length 150 \
  --seq-type paired \
  --meth-dir single_gene \
  --meth-spec one__Potri.001G055900.5.v3.0 \
  --circle-meth-layout separate \
  --circle-frame \
  --out Potri.001G055900.5.v3.0_circle_with_meth_separate.pdf
```

---

### 7. Parameter reference (grouped)
Required:
- `--gff3 PATH`, `--gene STR`, `--distance INT`, `--out PATH`

GFF3 key adapters:
- `--gene-keywords STR` (default: `ID`)
- `--exon-keywords STR` (default: `Parent`)

BAM rings:
- `--bam-dir DIR`
- `--bam-spec STR`
- `--name-spec STR` (must match)
- `--reads-length INT`
- `--seq-type {paired,single}`
- `--circle-bam-step INT`
- `--bam-group-ylims STR`

Methylation:
- `--meth-dir DIR`
- `--meth-spec STR`
- `--circle-meth-layout {combined,separate}`
- `--meth-drop-me-zero`

Geometry / consistency:
- `--circle-ring-start`, `--circle-ring-width`, `--circle-ring-gap`, `--circle-margin`, `--circle-lim`

Frame:
- `--circle-frame`, `--circle-frame-lw`

Figure size:
- `--fig-x`, `--fig-y`

---

### 8. Troubleshooting (EN)
- Gene not found: confirm transcript (mRNA) ID; adjust `--gene-keywords/--exon-keywords`.
- BAM slow/empty: ensure sorted + indexed BAM.
- Incomparable ring scales: use `--bam-group-ylims`.
- Inconsistent borders: set `--circle-lim` and keep geometry fixed.

