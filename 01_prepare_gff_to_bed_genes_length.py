#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""1trans_gff_to_bed_genes_length_EN.py

GFF3 -> BED6 + (CDS/exon) length table (NO --level)

What this script does
- Export gene / mRNA (or transcript) coordinates from a GFF3 file into BED6:
    chrom, start, end, ID, score, strand
- Summarize CDS or exon lengths.
  By default, lengths are first grouped by the attribute "Parent" (which usually points to transcript/mRNA).
- If BED IDs and length IDs mismatch AND --bed-feature=gene, the script will automatically try:
  transcript -> gene mapping via mRNA/transcript records (mRNA/transcript: ID=<tx>, Parent=<gene>)
  and then aggregate lengths to gene level.
- Supports multi-parent attributes (comma-separated) and will explode them.
- Optionally converts GFF3 coordinates (1-based inclusive) to standard BED coordinates (0-based half-open).

Typical use
- Build BED from mRNA + compute CDS length by transcript:
    --bed-feature mRNA --length-feature CDS
- Build BED from gene + compute CDS length by gene (CDS Parent is transcript):
    --bed-feature gene --length-feature CDS  (auto transcript->gene aggregation)

Notes
- GFF3 coordinates are 1-based inclusive; BED is 0-based half-open.
  If you enable --bed-zero-based, this script subtracts 1 from start, and keeps end unchanged.
"""

import os
import re
import sys
import argparse
import pandas as pd


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def extract_attr(series: pd.Series, key: str) -> pd.Series:
    """Extract attribute value from GFF3 attribute column by key (e.g., ID=xxx)."""
    pattern = rf"{re.escape(key)}=([^;]+)"
    return series.astype(str).str.extract(pattern, expand=False)


def explode_multi_value(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Split comma-separated values and explode into multiple rows."""
    df = df.copy()
    df[col] = df[col].fillna("").astype(str)
    df[col] = df[col].apply(lambda s: [x for x in s.split(",") if x] if s else [])
    df = df.explode(col)
    df = df[df[col].notna() & (df[col].astype(str) != "")]
    return df


def merge_intervals_sum(starts, ends) -> int:
    """Merge overlapping intervals and return total covered length (1-based inclusive)."""
    if len(starts) == 0:
        return 0
    pairs = sorted(zip(starts, ends), key=lambda x: (x[0], x[1]))
    cur_s, cur_e = pairs[0]
    total = 0
    for s, e in pairs[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += (cur_e - cur_s + 1)
            cur_s, cur_e = s, e
    total += (cur_e - cur_s + 1)
    return total


def summarize_length(df_feat: pd.DataFrame,
                     id_col: str,
                     start_col: str,
                     end_col: str,
                     merge_overlaps: bool) -> pd.DataFrame:
    """Summarize total length per id_col."""
    if df_feat.empty:
        return pd.DataFrame(columns=["length"]).astype({"length": int})

    if not merge_overlaps:
        tmp = df_feat.copy()
        tmp["length"] = tmp[end_col] - tmp[start_col] + 1
        return tmp.groupby(id_col, as_index=True)["length"].sum().to_frame()

    records = []
    for gid, sub in df_feat.groupby(id_col, sort=False):
        starts = sub[start_col].astype(int).tolist()
        ends = sub[end_col].astype(int).tolist()
        records.append((gid, merge_intervals_sum(starts, ends)))
    return pd.DataFrame(records, columns=[id_col, "length"]).set_index(id_col)


def intersection_size(a, b) -> int:
    return len(set(a).intersection(set(b)))


def main() -> None:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Generate a BED6 file and a CDS/exon length table from a GFF3 file. "
            "Auto-resolves ID mismatches (no --level)."
        )
    )

    # I/O
    parser.add_argument("-i", "--gff3", required=True,
                        help="Input GFF3 file path.")
    parser.add_argument("-o", "--out-bed", default="gene.bed",
                        help="Output BED6 file path.")
    parser.add_argument("-l", "--out-length", default="feature_length.tsv",
                        help="Output length table path (2 columns: ID<TAB>length).")

    # Feature selection
    parser.add_argument("--bed-feature", default="mRNA", choices=["gene", "mRNA", "transcript"],
                        help="Feature type used to build BED (gene, mRNA, or transcript).")
    parser.add_argument("--length-feature", default="CDS", choices=["CDS", "exon"],
                        help="Feature type used to compute length (CDS or exon).")

    # Attribute keys
    parser.add_argument("--id-key", default="ID",
                        help="Attribute key for feature ID (e.g., ID=xxx).")
    parser.add_argument("--parent-key", default="Parent",
                        help="Attribute key for parent ID (e.g., Parent=xxx).")

    # BED coordinate conversion
    parser.add_argument("--bed-zero-based", action="store_true",
                        help=(
                            "Convert GFF3 start (1-based inclusive) to BED start (0-based). "
                            "End stays unchanged (BED is half-open)."
                        ))

    # Other options
    parser.add_argument("--score", type=int, default=10,
                        help="Constant value for BED score column (5th column).")
    parser.add_argument("--merge-overlaps", action="store_true",
                        help="Merge overlapping CDS/exon segments before summing length (more accurate).")
    parser.add_argument("--min-match", type=int, default=10000,
                        help=(
                            "Minimum number of shared IDs between BED IDs and length IDs; "
                            "if lower, the script aborts to prevent wrong mapping. "
                            "(Lower this value for small datasets.)"
                        ))

    args = parser.parse_args()

    if not os.path.exists(args.gff3):
        raise FileNotFoundError(f"Input GFF3 not found: {args.gff3}")

    log(f"[INFO] Input GFF3: {args.gff3}")
    log(f"[INFO] BED feature: {args.bed_feature} | Length feature: {args.length_feature}")

    # Read GFF3 (only required columns)
    gff = pd.read_csv(
        args.gff3, sep="\t", comment="#", header=None,
        usecols=[0, 2, 3, 4, 6, 8],
        names=["ch", "ty", "st", "en", "zf", "attr"],
        dtype={"ch": str, "ty": str, "st": int, "en": int, "zf": str, "attr": str},
        low_memory=False
    )

    feature_types = sorted(gff["ty"].dropna().unique().tolist())
    log(f"[INFO] Feature types detected: {feature_types}")

    if args.bed_feature not in feature_types:
        raise ValueError(f"{args.bed_feature} not found in GFF3 feature types: {feature_types}")
    if args.length_feature not in feature_types:
        raise ValueError(f"{args.length_feature} not found in GFF3 feature types: {feature_types}")

    # 1) Build BED
    bed_df = gff[gff["ty"] == args.bed_feature].copy()
    bed_df["ID"] = extract_attr(bed_df["attr"], args.id_key)
    bed_df = bed_df[bed_df["ID"].notna()].copy()

    if args.bed_zero_based:
        bed_df["st"] = bed_df["st"] - 1

    bed_df["score"] = args.score
    bed = bed_df[["ch", "st", "en", "ID", "score", "zf"]].copy()

    bed_ids = bed["ID"].astype(str).tolist()
    log(f"[INFO] BED records: {len(bed)} | unique IDs: {len(set(bed_ids))}")

    # 2) Compute lengths by Parent
    feat = gff[gff["ty"] == args.length_feature].copy()
    feat["Parent"] = extract_attr(feat["attr"], args.parent_key)
    feat = feat[feat["Parent"].notna()].copy()
    feat = explode_multi_value(feat, "Parent")

    length_by_parent = summarize_length(
        df_feat=feat,
        id_col="Parent",
        start_col="st",
        end_col="en",
        merge_overlaps=args.merge_overlaps
    )
    length_by_parent.index = length_by_parent.index.astype(str)
    length_by_parent.index.name = "ID"

    parent_ids = length_by_parent.index.tolist()
    inter1 = intersection_size(bed_ids, parent_ids)
    log(f"[INFO] Intersection (BED IDs vs Parent IDs): {inter1}")

    length_out = None

    # Case 1: IDs already match
    if inter1 >= args.min_match:
        length_out = length_by_parent
        log("[INFO] Using length table grouped by Parent (IDs match BED).")

    # Case 2: mismatch and BED is gene -> try transcript->gene mapping
    elif args.bed_feature == "gene":
        tx_feature = None
        if "mRNA" in feature_types:
            tx_feature = "mRNA"
        elif "transcript" in feature_types:
            tx_feature = "transcript"

        if tx_feature is None:
            raise RuntimeError(
                "ID mismatch and cannot build transcript->gene mapping because neither 'mRNA' nor 'transcript' "
                "features are present in this GFF3.\n"
                "Suggestion: verify your GFF3 or check whether CDS/exon Parent refers directly to gene IDs."
            )

        tx = gff[gff["ty"] == tx_feature].copy()
        tx["tx_id"] = extract_attr(tx["attr"], args.id_key)          # transcript/mRNA ID
        tx["gene_id"] = extract_attr(tx["attr"], args.parent_key)    # Parent=gene ID
        tx = tx[tx["tx_id"].notna() & tx["gene_id"].notna()].copy()
        tx = explode_multi_value(tx, "gene_id")

        tx2gene = tx.set_index("tx_id")["gene_id"].astype(str)

        length_gene = length_by_parent.copy()
        length_gene["gene_id"] = length_gene.index.map(tx2gene)
        length_gene = length_gene[length_gene["gene_id"].notna()].copy()
        length_gene = length_gene.groupby("gene_id", as_index=True)["length"].sum().to_frame()
        length_gene.index = length_gene.index.astype(str)
        length_gene.index.name = "ID"

        inter2 = intersection_size(bed_ids, length_gene.index.tolist())
        log(f"[INFO] Intersection after transcript->gene aggregation: {inter2}")

        if inter2 >= args.min_match:
            length_out = length_gene
            log("[INFO] Using gene-aggregated length table (auto transcript->gene mapping).")
        else:
            raise RuntimeError(
                "ID mismatch persists after transcript->gene aggregation.\n"
                f"Intersection after aggregation ({inter2}) < --min-match ({args.min_match}).\n"
                "Possible reasons:\n"
                "- wrong --id-key / --parent-key\n"
                "- GFF3 uses different attribute keys (e.g., transcript_id/gene_id)\n"
                "- CDS/exon Parent refers to different IDs than expected\n"
            )

    # Case 3: mismatch for mRNA BED
    else:
        raise RuntimeError(
            "ID mismatch: BED is built from transcript/mRNA IDs, but CDS/exon Parent IDs do not match them.\n"
            f"Intersection ({inter1}) < --min-match ({args.min_match}).\n"
            "Suggestions:\n"
            "- Check whether CDS/exon Parent points to transcript/mRNA IDs\n"
            "- If Parent points to gene IDs, try --bed-feature gene\n"
            "- Verify --id-key / --parent-key for your GFF3"
        )

    # Write outputs
    bed.to_csv(args.out_bed, sep="\t", index=False, header=False)
    length_out.to_csv(args.out_length, sep="\t", index=True, header=False)

    log(f"[DONE] BED saved: {args.out_bed}")
    log(f"[DONE] Length table saved: {args.out_length}")


if __name__ == "__main__":
    main()
