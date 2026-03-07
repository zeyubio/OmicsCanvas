# 19_omicscanvas_go_enrich.py — GO Enrichment + Bubble Plot

This script performs **GO enrichment** using a **hypergeometric right‑tail test**, applies **Benjamini–Hochberg FDR**, and produces a **bubble scatter plot** (seaborn).

- Statistics: `p = P(X >= k),  X ~ Hypergeom(N, K, n)` fileciteturn49file0L22-L31
- Multiple testing: BH‑FDR (q‑value) fileciteturn49file0L30-L31

---

## 1) What it does

Given a query gene list and a GO annotation table:

1. Builds a GO‑annotated background (either **all annotated genes** or a user‑provided background list).
2. For each GO term:
   - `N`: number of background genes
   - `n`: number of query genes intersecting the background
   - `K`: number of background genes annotated with the term
   - `k`: number of query genes annotated with the term fileciteturn49file0L22-L28
3. Computes hypergeometric **right‑tail** p‑values (`sf(k-1, N, K, n)`).
4. Adjusts p‑values to **FDR** using BH. fileciteturn49file3L43-L67
5. Writes enrichment tables and draws bubble plots (Top‑N terms). fileciteturn49file1L1-L16

---

## 2) Requirements

Python packages:

- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `seaborn` fileciteturn49file0L58-L63

Install (pip):

```bash
pip install numpy pandas scipy matplotlib seaborn
```

---

## 3) Inputs

### 3.1 Query gene list (`--genes`) — required

Plain text, **one gene ID per line**. Empty lines and `#` comment lines are ignored. fileciteturn49file0L74-L90

Example:

```
GeneA
GeneB
GeneC
```

The script de‑duplicates genes while preserving the original order. fileciteturn49file0L84-L90

### 3.2 GO annotation table (`--go-annot`) — required

Tab‑separated text with **at least 4 columns**:

```
gene_id<TAB>GO_ID<TAB>Description<TAB>Ontology
```

Ontology must be one of:
- `biological_process`
- `molecular_function`
- `cellular_component` fileciteturn49file0L10-L13

Notes on parsing:
- Extra columns are allowed; only the first 4 columns are used. fileciteturn49file2L82-L84
- Header row is optional; the script tries to detect and drop it automatically. fileciteturn49file2L85-L87
- Lines starting with `#` are treated as comments.

### 3.3 Optional background gene list (`--background`)

If provided, background genes are `background_list ∩ (GO‑annotated genes)`. fileciteturn49file3L6-L10  
If omitted, background is **all genes present in GO annotations**. fileciteturn49file3L6-L8

Format: same as `--genes` (one gene per line).

---

## 4) Default filters

By default, the script **removes**:

- GO “root” terms: `GO:0008150`, `GO:0003674`, `GO:0005575` fileciteturn49file0L66-L71
- Terms with Description containing the word `obsolete` (case‑insensitive). fileciteturn49file3L1-L2

Override:
- `--keep-root` keeps root terms. fileciteturn49file4L48-L49
- `--keep-obsolete` keeps “obsolete” terms. fileciteturn49file4L49-L50

---

## 5) Outputs

### 5.1 Enrichment table (TSV)

For each run, the script writes:

- `<out_prefix>.<tag>.go_enrichment.tsv` fileciteturn49file1L1-L3 fileciteturn49file1L34-L36

Where `<tag>` is:
- `bp`, `cc`, `mf` when `--three-ontologies` fileciteturn49file1L18-L22
- otherwise: `all` or the selected ontology name (`biological_process`, …) fileciteturn49file1L23-L25

Columns include (core):
- `GO_ID`, `Description`, `Ontology`
- `k`, `n`, `K`, `N`, `GeneRatio`, `BgRatio`
- `pvalue`, `FDR`
- `geneIDs` (query genes hitting the term; `;`‑joined) fileciteturn49file3L45-L57

### 5.2 Bubble plot

If the enrichment result is non‑empty, the script also writes:

- `<out_prefix>.<tag>.go_bubbleplot.<pdf|png>` fileciteturn49file1L6-L14

Bubble plot encoding:
- **Y**: GO term Description (top terms)
- **X**: significance (`-log10(FDR)` or `-log10(pvalue)` by default)
- **Point size**: `GeneRatio`
- **Color**: `Ontology` fileciteturn49file4L9-L18

Legend is placed outside the axes; output is saved at `dpi=300`. fileciteturn49file4L24-L29

If no enriched terms pass filters, the plot is **not** generated. fileciteturn49file1L15-L16

---

## 6) Usage examples

### 6.1 Run BP/CC/MF separately (3 tables + 3 plots)

```bash
python 19_omicscanvas_go_enrich.py   --genes fig1_cluster_1_genes.txt   --go-annot test_GO.txt   --out-prefix cluster1   --three-ontologies   --top 30 --metric fdr --plot-format pdf
```
(Example mirrors the script header.) fileciteturn49file0L35-L42

### 6.2 Run a single ontology only

```bash
python 19_omicscanvas_go_enrich.py   --genes fig1_cluster_1_genes.txt   --go-annot test_GO.txt   --out-prefix cluster1   --ontology biological_process   --top 30 --metric fdr
```
(Example mirrors the script header.) fileciteturn49file0L43-L49

### 6.3 Use a custom background gene set

```bash
python 19_omicscanvas_go_enrich.py   --genes DEG_up.txt   --background all_detected_genes.txt   --go-annot go_annot.tsv   --out-prefix DEG_up.bg_detected   --ontology all   --top 30 --metric fdr
```

### 6.4 Plot by pvalue and do not transform X

```bash
python 19_omicscanvas_go_enrich.py   --genes DEG_up.txt   --go-annot go_annot.tsv   --out-prefix DEG_up.pvalue_raw   --metric pvalue   --x-transform none   --plot-format png
```

### 6.5 Keep root terms / obsolete terms

```bash
python 19_omicscanvas_go_enrich.py   --genes DEG_up.txt   --go-annot go_annot.tsv   --out-prefix DEG_up.keep_root   --keep-root --keep-obsolete   --ontology biological_process
```

---

## 7) Arguments

Required:
- `--genes` : query gene list file fileciteturn49file4L38-L39
- `--go-annot` : GO annotation TSV (>=4 columns) fileciteturn49file4L39-L40
- `--out-prefix` : output prefix fileciteturn49file4L50-L50

Ontology selection:
- `--three-ontologies` : run BP/CC/MF separately (ignores `--ontology`) fileciteturn49file4L55-L59
- `--ontology {all, biological_process, molecular_function, cellular_component}` (default `all`) fileciteturn49file4L40-L45

Background and filters:
- `--background` : optional background gene list (default: all GO‑annotated genes) fileciteturn49file4L46-L47
- `--min-k` (default `2`) : minimum overlap genes (`k`) for a term to be kept fileciteturn49file4L47-L48
- `--keep-root` : keep root terms fileciteturn49file4L48-L49
- `--keep-obsolete` : keep “obsolete” terms fileciteturn49file4L49-L50

Plot settings:
- `--top` (default `30`) : top terms to show fileciteturn49file4L51-L51
- `--metric {pvalue,fdr}` (default `fdr`) : x-axis significance metric fileciteturn49file4L52-L53
- `--x-transform {neglog10,none}` (default `neglog10`) : x-axis transform fileciteturn49file4L53-L54
- `--plot-format {pdf,png}` (default `pdf`) : plot output format fileciteturn49file4L54-L55

---

## 8) Troubleshooting

### “None of the query genes are in the GO-annotated background”
Your query gene IDs do not overlap with GO‑annotated gene IDs (or with your `--background` after intersecting). fileciteturn49file3L15-L18  
Check:
- gene ID versioning (e.g., `.1` suffix vs not)
- whether your GO annotation uses the same ID system as your gene list

### “Background gene set is empty…”
Your `--background` has no overlap with GO‑annotated genes. fileciteturn49file3L11-L14  
Fix by using a compatible background list, or omit `--background`.

### No plot generated
If all terms are filtered out (e.g., `--min-k` too high, or few genes), the enrichment table may be empty and the script will skip plotting. fileciteturn49file1L15-L16  
Try:
- decreasing `--min-k`
- increasing `--top`
- ensuring enough query genes are annotated

---

## 9) CLI help

```bash
python 19_omicscanvas_go_enrich.py -h
```
