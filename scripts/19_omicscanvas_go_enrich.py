#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OmicsCanvas - GO enrichment (hypergeometric) + BH-FDR + bubble scatter plot (seaborn)

Inputs
------
1) --genes: Gene list file, one gene ID per line
2) --go-annot: GO annotation TSV with at least 4 columns:
      gene_id<TAB>GO_ID<TAB>Description<TAB>Ontology
   Ontology should be one of:
      biological_process, molecular_function, cellular_component

Outputs
-------
- <out_prefix>.<bp|cc|mf>.go_enrichment.tsv  (when --three-ontologies)
- <out_prefix>.<bp|cc|mf>.go_bubbleplot.(pdf|png)

Notes
-----
Statistics:
  Hypergeometric right-tail p-value:
    N = number of background genes
    n = number of query genes (intersect background)
    K = number of background genes annotated with a term
    k = number of query genes annotated with a term
    p = P(X >= k), X ~ Hypergeom(N, K, n)

Multiple testing:
  Benjamini–Hochberg FDR (BH-qvalue)

Examples
--------
Run BP/CC/MF separately (Top30 each):
  python omicscanvas_go_enrich.py \
    --genes fig1_cluster_1_genes.txt \
    --go-annot test_GO.txt \
    --out-prefix cluster1 \
    --three-ontologies \
    --top 30 --metric fdr --plot-format pdf

Run a single ontology only:
  python omicscanvas_go_enrich.py \
    --genes fig1_cluster_1_genes.txt \
    --go-annot test_GO.txt \
    --out-prefix cluster1 \
    --ontology biological_process \
    --top 30 --metric fdr
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

import matplotlib.pyplot as plt
import seaborn as sns


ROOT_TERMS = {
    "GO:0008150",  # biological_process
    "GO:0003674",  # molecular_function
    "GO:0005575",  # cellular_component
}
ROOT_NAMES = {"biological_process", "molecular_function", "cellular_component"}


def read_gene_list(path: str) -> list[str]:
    """Read one-gene-per-line text, return unique genes while preserving order."""
    genes: list[str] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            genes.append(s)

    seen = set()
    uniq = []
    for g in genes:
        if g not in seen:
            uniq.append(g)
            seen.add(g)
    return uniq


def _looks_like_header(row0: pd.Series) -> bool:
    """Heuristic: detect header-like first row."""
    vals = [str(x).strip().lower() for x in row0.tolist()]
    joined = " ".join(vals)
    if "go" in joined and "gene" in joined:
        return True
    if vals and (vals[0] in {"gene", "gene_id", "id"}):
        return True
    return False


def parse_go_annot(path: str) -> pd.DataFrame:
    """
    Parse GO annotation TSV (>=4 columns):
      gene, go_id, description, ontology
    Robust to extra columns; ignores malformed lines.
    """
    # Try fast path with pandas
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            comment="#",
            dtype=str,
            keep_default_na=False,
            na_values=[],
        )
        if df.shape[1] < 4:
            raise ValueError("GO annotation must have >=4 columns (TSV).")

        df = df.iloc[:, :4].copy()
        df.columns = ["gene", "go_id", "description", "ontology"]

        if df.shape[0] > 0 and _looks_like_header(df.iloc[0]):
            df = df.iloc[1:].copy()

    except Exception:
        # Fallback manual parsing
        rows = []
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                rows.append((parts[0], parts[1], parts[2], parts[3]))
        df = pd.DataFrame(rows, columns=["gene", "go_id", "description", "ontology"])

    # cleanup
    df["gene"] = df["gene"].astype(str).str.strip()
    df["go_id"] = df["go_id"].astype(str).str.strip()
    df["description"] = df["description"].astype(str).str.strip()
    df["ontology"] = df["ontology"].astype(str).str.strip()

    df = df[(df["gene"] != "") & (df["go_id"] != "") & (df["description"] != "") & (df["ontology"] != "")]
    return df.reset_index(drop=True)


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR adjustment.
    Returns q-values aligned to original order.
    """
    pvals = np.asarray(pvals, dtype=float)
    m = pvals.size
    if m == 0:
        return pvals

    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * m / (np.arange(1, m + 1))

    # enforce monotonicity from end
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    out = np.empty_like(q)
    out[order] = q
    return out


def go_enrichment(
    query_genes: list[str],
    go_df: pd.DataFrame,
    ontology: str = "all",
    background_genes: list[str] | None = None,
    drop_root: bool = True,
    drop_obsolete: bool = True,
    min_k: int = 2,
) -> pd.DataFrame:
    """
    Hypergeometric enrichment + BH-FDR.
    """
    df = go_df.copy()

    if ontology != "all":
        df = df[df["ontology"] == ontology].copy()

    if drop_root:
        df = df[~df["go_id"].isin(ROOT_TERMS)].copy()
        df = df[~df["description"].isin(ROOT_NAMES)].copy()

    if drop_obsolete:
        df = df[~df["description"].str.contains(r"\bobsolete\b", case=False, na=False)].copy()

    annotated_genes = set(df["gene"].unique().tolist())

    if background_genes is None:
        bg_genes = annotated_genes
    else:
        bg_genes = set(background_genes) & annotated_genes

    N = len(bg_genes)
    if N == 0:
        raise ValueError("Background gene set is empty after intersecting with GO-annotated genes.")

    query_set = set(query_genes) & bg_genes
    n = len(query_set)
    if n == 0:
        raise ValueError("None of the query genes are in the GO-annotated background.")

    df_bg = df[df["gene"].isin(bg_genes)].copy()

    # term -> background genes
    term2genes_bg: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for r in df_bg.itertuples(index=False):
        term2genes_bg[(r.go_id, r.description, r.ontology)].add(r.gene)

    # term -> query genes
    df_q = df_bg[df_bg["gene"].isin(query_set)]
    term2genes_q: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for r in df_q.itertuples(index=False):
        term2genes_q[(r.go_id, r.description, r.ontology)].add(r.gene)

    rows = []
    for key, bg_term_genes in term2genes_bg.items():
        go_id, desc, ont = key
        K = len(bg_term_genes)
        q_term_genes = term2genes_q.get(key, set())
        k = len(q_term_genes)

        if k < min_k:
            continue

        pval = float(hypergeom.sf(k - 1, N, K, n))
        rows.append(
            {
                "GO_ID": go_id,
                "Description": desc,
                "Ontology": ont,
                "k": k,
                "n": n,
                "K": K,
                "N": N,
                "GeneRatio": k / n,
                "BgRatio": K / N,
                "pvalue": pval,
                "geneIDs": ";".join(sorted(q_term_genes)),
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=["GO_ID", "Description", "Ontology", "k", "n", "K", "N", "GeneRatio", "BgRatio", "pvalue", "FDR", "geneIDs"]
        )

    res = pd.DataFrame(rows).sort_values("pvalue", ascending=True).reset_index(drop=True)
    res["FDR"] = benjamini_hochberg(res["pvalue"].values)
    return res


def bubble_plot(
    res: pd.DataFrame,
    outpath: str,
    metric: str = "fdr",
    top: int = 30,
    x_transform: str = "neglog10",
    figsize_w: float = 9.0,
):
    if res.empty:
        raise ValueError("No enrichment results to plot (empty table).")

    metric = metric.lower()
    if metric not in {"pvalue", "fdr"}:
        raise ValueError("--metric must be 'pvalue' or 'fdr'.")

    xcol = "pvalue" if metric == "pvalue" else "FDR"

    plot_df = res.sort_values(xcol, ascending=True).head(int(top)).copy()

    if x_transform == "none":
        plot_df["_x"] = plot_df[xcol].astype(float)
        xlab = xcol
    elif x_transform == "neglog10":
        plot_df["_x"] = -np.log10(np.clip(plot_df[xcol].astype(float), 1e-300, 1.0))
        xlab = f"-log10({xcol})"
    else:
        raise ValueError("--x-transform must be 'neglog10' or 'none'.")

    # order y by significance (top at top)
    plot_df = plot_df.sort_values("_x", ascending=False)
    plot_df["Description"] = plot_df["Description"].astype(str)
    plot_df["Description"] = pd.Categorical(plot_df["Description"], categories=plot_df["Description"].tolist(), ordered=True)

    height = max(3.5, 0.35 * len(plot_df) + 1.2)
    plt.figure(figsize=(figsize_w, height))

    ax = sns.scatterplot(
        data=plot_df,
        x="_x",
        y="Description",
        size="GeneRatio",
        hue="Ontology",
        sizes=(40, 420),
        alpha=0.9,
        linewidth=0.3,
    )

    ax.set_xlabel(xlab)
    ax.set_ylabel("")
    ax.grid(True, axis="x", linestyle="--", alpha=0.25)

    # place legend outside
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0, frameon=True)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        prog="omicscanvas_go_enrich.py",
        description="GO enrichment (hypergeometric) + BH-FDR + bubble scatter plot (seaborn).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--genes", required=True, help="Gene list file (one gene ID per line).")
    ap.add_argument("--go-annot", required=True, help="GO annotation TSV: gene<TAB>GO_ID<TAB>Description<TAB>Ontology")
    ap.add_argument(
        "--ontology",
        default="all",
        choices=["all", "biological_process", "molecular_function", "cellular_component"],
        help="GO namespace to test (ignored if --three-ontologies is set).",
    )
    ap.add_argument("--background", default=None, help="Optional background gene list. Default: all GO-annotated genes.")
    ap.add_argument("--min-k", type=int, default=2, help="Minimum overlap genes (k) to keep a GO term.")
    ap.add_argument("--keep-root", action="store_true", help="Keep root terms (GO:0008150/GO:0003674/GO:0005575).")
    ap.add_argument("--keep-obsolete", action="store_true", help="Keep terms whose Description contains 'obsolete'.")
    ap.add_argument("--out-prefix", required=True, help="Output prefix.")
    ap.add_argument("--top", type=int, default=30, help="Top terms to show in bubble plot.")
    ap.add_argument("--metric", choices=["pvalue", "fdr"], default="fdr", help="X-axis significance metric.")
    ap.add_argument("--x-transform", choices=["neglog10", "none"], default="neglog10", help="Transform x for plotting.")
    ap.add_argument("--plot-format", choices=["pdf", "png"], default="pdf", help="Bubble plot format.")
    ap.add_argument(
        "--three-ontologies",
        action="store_true",
        help="Generate BP/CC/MF separately (3 tables + 3 plots).",
    )
    args = ap.parse_args()

    query_genes = read_gene_list(args.genes)
    go_df = parse_go_annot(args.go_annot)

    bg_genes = None
    if args.background is not None:
        bg_genes = read_gene_list(args.background)

    n_query = len(set(query_genes))
    n_annot_query = len(set(query_genes) & set(go_df["gene"].unique()))
    print(f"[INFO] Query genes: {n_query}; annotated-in-background: {n_annot_query}")

    def run_one(ont: str, tag: str):
        res = go_enrichment(
            query_genes=query_genes,
            go_df=go_df,
            ontology=ont,
            background_genes=bg_genes,
            drop_root=not args.keep_root,
            drop_obsolete=not args.keep_obsolete,
            min_k=args.min_k,
        )

        out_tsv = f"{args.out_prefix}.{tag}.go_enrichment.tsv"
        res.to_csv(out_tsv, sep="\t", index=False)
        print(f"[INFO] Wrote: {out_tsv}")

        if not res.empty:
            out_plot = f"{args.out_prefix}.{tag}.go_bubbleplot.{args.plot_format}"
            bubble_plot(
                res=res,
                outpath=out_plot,
                metric=args.metric,
                top=args.top,
                x_transform=args.x_transform,
            )
            print(f"[INFO] Wrote: {out_plot}")
        else:
            print(f"[WARN] No enriched terms for {tag}; plot not generated.")

    if args.three_ontologies:
        run_one("biological_process", "bp")
        run_one("cellular_component", "cc")
        run_one("molecular_function", "mf")
    else:
        # single ontology / all
        tag = args.ontology if args.ontology != "all" else "all"
        res = go_enrichment(
            query_genes=query_genes,
            go_df=go_df,
            ontology=args.ontology,
            background_genes=bg_genes,
            drop_root=not args.keep_root,
            drop_obsolete=not args.keep_obsolete,
            min_k=args.min_k,
        )
        out_tsv = f"{args.out_prefix}.{tag}.go_enrichment.tsv"
        res.to_csv(out_tsv, sep="\t", index=False)
        print(f"[INFO] Wrote: {out_tsv}")

        if not res.empty:
            out_plot = f"{args.out_prefix}.{tag}.go_bubbleplot.{args.plot_format}"
            bubble_plot(
                res=res,
                outpath=out_plot,
                metric=args.metric,
                top=args.top,
                x_transform=args.x_transform,
            )
            print(f"[INFO] Wrote: {out_plot}")
        else:
            print("[WARN] No enriched terms passed filters; plot not generated.")


if __name__ == "__main__":
    main()
