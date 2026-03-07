#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""BamTrackSuite: CX context splitter

Purpose
-------
Split a Bismark CX_report into three context-specific coverage files (CG / CHG / CHH).
Each output file contains four columns:
  chrom  position  methylated_count  total_count

Rationale
---------
This script is intentionally minimal: it performs context splitting and optional depth
capping only. The output is designed to be consumed by downstream BAM/track-centric
visualization or multi-omics integration workflows.

Input
-----
Bismark CX_report (tab-delimited, typically no header). A common column layout is:
  1) chromosome
  2) position (1-based)
  3) strand
  4) count_methylated
  5) count_unmethylated
  6) context (CG/CHG/CHH/...)

This script uses: chromosome, position, count_methylated, count_unmethylated, context.

Output
------
Three files are written into the output directory:
  <prefix>_CG.CX
  <prefix>_CHG.CX
  <prefix>_CHH.CX

Arguments (args)
----------------
-i/--input-cx : str
    Input CX_report path.
-p/--prefix : str
    Output file prefix (usually a sample name).
-d/--out-dir : str
    Output directory (created if missing).
--max-depth : int
    Depth cap for total_count. If total_count exceeds max_depth, it is truncated.
    methylated_count is additionally constrained to never exceed total_count.
--chunksize : int
    Rows per chunk for streaming large CX_report files.

Example
-------
python cx_context_split.en.py \
  -i sample.CX_report.txt \
  -p sample \
  -d meth_data \
  --max-depth 300 \
  --chunksize 1000000
"""

from __future__ import annotations

import argparse
import os
from typing import Tuple

import pandas as pd


SUITE_NAME = "BamTrackSuite"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bts-cx-split",
        description=f"{SUITE_NAME} | Split Bismark CX_report into CG/CHG/CHH coverage files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input-cx",
        required=True,
        metavar="FILE",
        help="Input Bismark CX_report file path (TSV, no header).",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        required=True,
        metavar="STR",
        help=(
            "Output file prefix (sample name). Outputs: <prefix>_CG.CX, <prefix>_CHG.CX, <prefix>_CHH.CX."
        ),
    )
    parser.add_argument(
        "-d",
        "--out-dir",
        default="meth_data",
        metavar="DIR",
        help="Output directory (created if missing).",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=300,
        metavar="INT",
        help=(
            "Depth cap: total_count > max_depth will be truncated to max_depth; methylated_count is constrained to <= total_count."
        ),
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=1_000_000,
        metavar="INT",
        help="Rows per chunk for streaming large CX_report files.",
    )

    parser.epilog = (
        "Example:\n"
        "  bts-cx-split -i sample.CX_report.txt -p sample -d meth_data --max-depth 300 --chunksize 1000000\n"
        "Outputs will be meth_data/sample_CG.CX, etc."
    )

    return parser


def ensure_output_paths(out_dir: str, prefix: str) -> Tuple[str, str, str]:
    os.makedirs(out_dir, exist_ok=True)

    cg_path = os.path.join(out_dir, f"{prefix}_CG.CX")
    chg_path = os.path.join(out_dir, f"{prefix}_CHG.CX")
    chh_path = os.path.join(out_dir, f"{prefix}_CHH.CX")

    # Remove existing outputs to avoid accidental appends
    for p in (cg_path, chg_path, chh_path):
        if os.path.exists(p):
            os.remove(p)

    return cg_path, chg_path, chh_path


def split_cx_report(
    input_cx: str,
    prefix: str,
    out_dir: str,
    max_depth: int,
    chunksize: int,
) -> Tuple[str, str, str]:
    if not os.path.isfile(input_cx):
        raise FileNotFoundError(f"Input file not found: {input_cx}")

    cg_path, chg_path, chh_path = ensure_output_paths(out_dir, prefix)

    # CX_report typically has no header; we load columns [0,1,3,4,5]
    reader = pd.read_csv(
        input_cx,
        sep="\t",
        usecols=[0, 1, 3, 4, 5],
        header=None,
        names=["chrom", "pos", "me", "un", "ctx"],
        chunksize=chunksize,
        dtype={"chrom": str, "pos": int, "me": int, "un": int, "ctx": str},
    )

    total_rows = 0
    for _, chunk in enumerate(reader, start=1):
        total_rows += len(chunk)

        chunk["total"] = (chunk["me"] + chunk["un"]).astype(int)

        meth = chunk.loc[chunk["total"] != 0, ["chrom", "pos", "me", "total", "ctx"]].copy()

        # Cap depth (cap total first, then enforce me<=total)
        if max_depth is not None and max_depth > 0:
            meth.loc[meth["total"] > max_depth, "total"] = max_depth
            meth.loc[meth["me"] > meth["total"], "me"] = meth["total"]

        for ctx, out_path in (("CG", cg_path), ("CHG", chg_path), ("CHH", chh_path)):
            sub = meth.loc[meth["ctx"] == ctx, ["chrom", "pos", "me", "total"]]
            if not sub.empty:
                sub.to_csv(out_path, sep="\t", index=False, header=False, mode="a")

    print(f"[{SUITE_NAME}] Done. Total rows processed: {total_rows}")
    print(f"  CG : {cg_path}")
    print(f"  CHG: {chg_path}")
    print(f"  CHH: {chh_path}")

    return cg_path, chg_path, chh_path


def main() -> None:
    args = build_parser().parse_args()
    split_cx_report(
        input_cx=args.input_cx,
        prefix=args.prefix,
        out_dir=args.out_dir,
        max_depth=args.max_depth,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()
