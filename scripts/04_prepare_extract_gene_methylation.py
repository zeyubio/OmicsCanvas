#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""omicscanvas_extract_gene_methylation_EN.py

Extract per-base methylation (CX) records for one gene (or a gene list) from
precomputed cytosine-context files (CX).

This script is designed for downstream OmicsCanvas plotting/quantification.
It does NOT change your methylation calling; it only subsets CX records by
(gene coordinates ± distance) and converts genomic coordinates to a
strand-aware *relative position*.

Input CX file format (no header by default):
    ch\tpo\tme\tal
where:
    - ch: chromosome name (string)
    - po: 1-based genomic position (integer)
    - me: methylated counts (or methylation numerator)
    - al: total counts (or methylation denominator)

Relative coordinate convention:
    - For '+' strand: rel_po = po - (gene_start - distance)
    - For '-' strand: rel_po = (gene_end + distance) - po
So rel_po always increases in the transcription direction (5' -> 3').

Outputs (TSV):
    name\tpo\tme\tal
"""

import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Optional, List


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="omicscanvas_extract_gene_methylation.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Extract per-base methylation (CX) records for one gene or a gene list.\n\n"
            "Inputs:\n"
            "  1) GFF3 annotation (to get gene/transcript coordinates)\n"
            "  2) CX files under --cx-dir, named: <SRR>_<CX><cx-suffix>\n\n"
            "Outputs:\n"
            "  One TSV per (SRR/prefix, gene, context) in --outdir\n"
            "  Columns: name, po (relative position), me, al\n"
        ),
    )

    # --- Annotation ---
    p.add_argument(
        "-g",
        "--gff3",
        required=True,
        help=(
            "Input GFF3 file (tab-delimited). Example:\n"
            "  genome/Ptrichocarpa_210_v3.0.gene.gff3"
        ),
    )
    p.add_argument(
        "--feature-type",
        default="mRNA",
        help=(
            "Which GFF3 feature type to use as the target coordinate. (default: mRNA)\n"
            "Common choices: gene, mRNA, transcript.\n"
            "Tip: If your IDs are on gene lines, use --feature-type gene."
        ),
    )
    p.add_argument(
        "--attr-key",
        default="ID",
        help=(
            "Which attribute key to extract the identifier from the GFF3 attribute column. (default: ID)\n"
            "Example: ID, Parent, gene_id, transcript_id.\n"
            "The script will search 'attr-key=<value>' in the 9th column."
        ),
    )

    # --- CX files ---
    p.add_argument(
        "--cx-dir",
        required=True,
        help=(
            "Directory containing CX files. Example: ./meth/meth_data/\n"
            "Each file is expected as: <SRR>_<CX><cx-suffix>\n"
            "Example: SRR9321764_CHH.CX"
        ),
    )
    p.add_argument(
        "--cx-suffix",
        default=".CX",
        help=(
            "CX filename suffix (default: .CX).\n"
            "If your files are compressed, you may use .CX.gz (pandas can read gz)."
        ),
    )
    p.add_argument(
        "--contexts",
        default="CG,CHG,CHH",
        help=(
            "Contexts to extract, comma-separated (default: CG,CHG,CHH).\n"
            "Each context corresponds to one CX file: <SRR>_<context><cx-suffix>."
        ),
    )
    p.add_argument(
        "--srr",
        required=True,
        help=(
            "SRR list, comma-separated. Example: SRR9321764,SRRXXXX\n"
            "For each SRR, the script will look for <SRR>_<CX><cx-suffix>."
        ),
    )
    p.add_argument(
        "--prefix",
        default=None,
        help=(
            "Optional output prefix list (comma-separated).\n"
            "If not set, prefixes are identical to SRR IDs.\n"
            "Length can be shorter than SRR list; missing prefixes will be auto-filled by SRR."
        ),
    )

    # --- Gene selection ---
    p.add_argument(
        "--gene",
        default=None,
        help=(
            "Single target gene ID (exact match to the extracted GFF3 ID).\n"
            "Example: Potri.001G055900.5.v3.0"
        ),
    )
    p.add_argument(
        "--gene-list",
        default=None,
        help=(
            "Gene list file, one gene ID per line. Lines starting with '#' are ignored.\n"
            "If both --gene and --gene-list are provided, --gene-list is used."
        ),
    )

    # --- Window and IO ---
    p.add_argument(
        "--distance",
        type=int,
        default=1000,
        help=(
            "Upstream/downstream window size in bp. (default: 1000)\n"
            "Extract region: [start-distance, end+distance].\n"
            "Note: start-distance is clipped to >= 1."
        ),
    )
    p.add_argument(
        "-o",
        "--outdir",
        default="single_gene",
        help="Output directory (default: single_gene).",
    )
    p.add_argument(
        "--out-suffix",
        default=".tsv",
        help="Output file suffix (default: .tsv).",
    )
    p.add_argument(
        "--chunksize",
        type=int,
        default=2_000_000,
        help=(
            "Read CX file in chunks (rows) to reduce memory. (default: 2000000)\n"
            "Set 0 to read the entire file at once (faster but uses more RAM)."
        ),
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output files if they already exist.",
    )

    return p.parse_args()


def read_gene_ids(gene: Optional[str], gene_list: Optional[str]) -> List[str]:
    """Return a list of gene IDs to process."""
    if gene_list:
        genes: List[str] = []
        with open(gene_list, "r", encoding="utf-8") as f:
            for line in f:
                x = line.strip()
                if x and not x.startswith("#"):
                    genes.append(x)
        genes = list(dict.fromkeys(genes))  # unique, keep order
        if not genes:
            raise ValueError(f"Gene list is empty: {gene_list}")
        return genes

    if gene:
        return [gene.strip()]

    raise ValueError("You must provide --gene or --gene-list")


def safe_filename(x: str) -> str:
    """Convert any string into a safe filename token."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", x)


def load_gff_targets(
    gff3_path: str,
    feature_type: str,
    attr_key: str,
    genes: set,
    distance: int,
) -> pd.DataFrame:
    """Load target features from GFF3 and extend them by distance."""
    gff = pd.read_csv(
        gff3_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 2, 3, 4, 6, 8],
        names=["ch", "ty", "st", "en", "zf", "attr"],
        dtype={"ch": str, "ty": str, "st": np.int64, "en": np.int64, "zf": str, "attr": str},
        low_memory=False,
    )

    pat = rf"{re.escape(attr_key)}=([^;]+)"
    sub = gff[gff["ty"] == feature_type].copy()
    sub["ID"] = sub["attr"].str.extract(pat, expand=False)

    sub = sub[sub["ID"].isin(genes)].copy()
    if sub.empty:
        raise ValueError(
            f"No target IDs found in GFF3 (feature={feature_type}, attr_key={attr_key})."
        )

    sub["po1"] = (sub["st"] - distance).clip(lower=1)
    sub["po2"] = sub["en"] + distance

    sub = sub[["ch", "po1", "po2", "zf", "ID"]].drop_duplicates()
    sub.sort_values(["ch", "po1"], inplace=True)
    return sub


def read_cx_subset(cx_file: str, chrom: str, start: int, end: int, chunksize: int) -> pd.DataFrame:
    """Read a chromosome+interval subset from one CX file."""
    usecols = [0, 1, 2, 3]
    names = ["ch", "po", "me", "al"]
    dtypes = {"ch": str, "po": np.int64, "me": np.float64, "al": np.float64}

    if chunksize and chunksize > 0:
        it = pd.read_csv(
            cx_file,
            sep="\t",
            header=None,
            names=names,
            usecols=usecols,
            dtype=dtypes,
            chunksize=chunksize,
            low_memory=False,
        )
        keep = []
        for chunk in it:
            c = chunk[chunk["ch"] == chrom]
            if not c.empty:
                c = c[(c["po"] >= start) & (c["po"] <= end)]
                if not c.empty:
                    keep.append(c)
        if keep:
            out = pd.concat(keep, ignore_index=True)
        else:
            out = pd.DataFrame(columns=names)
    else:
        df = pd.read_csv(
            cx_file,
            sep="\t",
            header=None,
            names=names,
            usecols=usecols,
            dtype=dtypes,
            low_memory=False,
        )
        out = df[(df["ch"] == chrom) & (df["po"] >= start) & (df["po"] <= end)].copy()

    return out


def extract_one_gene(meth_df: pd.DataFrame, po1: int, po2: int, strand: str, gene_id: str) -> pd.DataFrame:
    """Convert genomic positions to strand-aware relative positions for one gene window."""
    if meth_df.empty:
        return pd.DataFrame(columns=["name", "po", "me", "al"])

    po = meth_df["po"].to_numpy(np.int64)
    me = meth_df["me"].to_numpy()
    al = meth_df["al"].to_numpy()

    mask = (po >= po1) & (po <= po2)
    if not np.any(mask):
        return pd.DataFrame(columns=["name", "po", "me", "al"])

    po = po[mask]
    me = me[mask]
    al = al[mask]

    if strand == "+":
        rel = po - po1
    else:
        rel = po2 - po

    out = pd.DataFrame({"name": gene_id, "po": rel.astype(np.int64), "me": me, "al": al})
    out.sort_values("po", inplace=True)
    return out


def main() -> None:
    args = parse_args()

    cx_dir = Path(args.cx_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = read_gene_ids(args.gene, args.gene_list)
    gene_set = set(genes)

    # SRR / prefix
    srr_list = [x.strip() for x in args.srr.split(",") if x.strip()]
    if not srr_list:
        raise ValueError("--srr is empty")

    if args.prefix:
        prefix_list = [x.strip() for x in args.prefix.split(",") if x.strip()]
        while len(prefix_list) < len(srr_list):
            prefix_list.append(srr_list[len(prefix_list)])
    else:
        prefix_list = srr_list[:]

    contexts = [x.strip() for x in args.contexts.split(",") if x.strip()]
    if not contexts:
        raise ValueError("--contexts is empty")

    # Load targets
    targets = load_gff_targets(args.gff3, args.feature_type, args.attr_key, gene_set, args.distance)

    found_genes = set(targets["ID"].unique().tolist())
    missing = [g for g in genes if g not in found_genes]
    if missing:
        print(f"[WARN] {len(missing)} gene(s) not found in GFF3: " + ", ".join(missing))

    # Group by chromosome
    by_chrom = {}
    for chrom, dfc in targets.groupby("ch", sort=False):
        start = int(dfc["po1"].min())
        end = int(dfc["po2"].max())
        by_chrom[chrom] = (dfc.reset_index(drop=True), start, end)

    # Main loop
    for srr, prefix in zip(srr_list, prefix_list):
        for cx in contexts:
            cx_file = cx_dir / f"{srr}_{cx}{args.cx_suffix}"
            if not cx_file.exists():
                print(f"[WARN] Missing file: {cx_file}")
                continue

            print(f"[INFO] Processing: SRR={srr}  CX={cx}  file={cx_file.name}")

            for chrom, (dfc, start, end) in by_chrom.items():
                meth = read_cx_subset(str(cx_file), chrom, start, end, args.chunksize)

                for _, row in dfc.iterrows():
                    gene_id = row["ID"]
                    po1 = int(row["po1"])
                    po2 = int(row["po2"])
                    zf = row["zf"]

                    out_df = extract_one_gene(meth, po1, po2, zf, gene_id)

                    out_name = f"{safe_filename(prefix)}__{safe_filename(gene_id)}__{cx}{args.out_suffix}"
                    out_path = outdir / out_name

                    if out_path.exists() and (not args.overwrite):
                        continue

                    out_df.to_csv(out_path, sep="\t", index=False)
                    print(f"[OK] Saved: {out_path}")

    print("[DONE]")


if __name__ == "__main__":
    main()
