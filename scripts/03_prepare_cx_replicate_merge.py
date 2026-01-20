#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""BamTrackSuite: CX replicate merger

Purpose
-------
Merge 2 or 3 replicate coverage files by (chrom, position) and sum column 3 and 4.
This is commonly used to merge CG/CHG/CHH coverage files produced from a Bismark
CX_report, but it also works for any file with the same first 4 columns:
  chrom  pos  col3  col4

Input requirements
------------------
- Each input must be sorted by (chrom, pos) in ascending order.
- Plain text or gzip-compressed text (.gz) is supported.
- Each line must have at least 4 columns; extra columns are ignored.

Output
------
A 4-column file:
  chrom  pos  sum_col3  sum_col4
If the output path ends with .gz, gzip output is written.

Arguments (args)
----------------
-o/--out : str
    Output file path. If it ends with .gz, gzip output is produced.
inputs : list[str]
    2 or 3 input replicate files (.tsv/.txt or .gz). The order does not matter,
    but each file must be internally sorted.

Examples
--------
python cx_replicate_merge.en.py -o merged_CG.CX rep1_CG.CX rep2_CG.CX
python cx_replicate_merge.en.py -o merged_CG.CX.gz rep1.gz rep2.gz rep3.gz
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
from typing import Optional, Tuple, List, TextIO


SUITE_NAME = "BamTrackSuite"

_chr_re = re.compile(r"^([A-Za-z_]+)(\d+)$")


def open_text_read(path: str) -> TextIO:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def open_text_write(path: str) -> TextIO:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "wb"), encoding="utf-8")
    return open(path, "w", encoding="utf-8")


def chr_sort_key(chrom: str) -> Tuple:
    """Chromosome ordering helper for key comparison (inputs are assumed sorted)."""
    m = _chr_re.match(chrom)
    if m:
        prefix, num = m.group(1), int(m.group(2))
        return (prefix, num)
    return (chrom,)


def parse_line(line: str) -> Optional[Tuple[str, int, int, int]]:
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    parts = line.split()
    if len(parts) < 4:
        raise ValueError(f"Bad line (need >=4 columns): {line}")
    chrom = parts[0]
    pos = int(parts[1])
    v3 = int(parts[2])
    v4 = int(parts[3])
    return chrom, pos, v3, v4


def read_next(handle: TextIO) -> Optional[Tuple[Tuple, str, int, int, int]]:
    for line in handle:
        parsed = parse_line(line)
        if parsed is None:
            continue
        chrom, pos, v3, v4 = parsed
        key = (chr_sort_key(chrom), pos)
        return key, chrom, pos, v3, v4
    return None


def min_key(entries: List[Optional[Tuple]]) -> Optional[Tuple]:
    keys = [e[0] for e in entries if e is not None]
    return min(keys) if keys else None


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="bts-cx-merge",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=f"{SUITE_NAME} | Merge 2/3 replicate files by (chrom,pos) and sum columns 3/4",
    )
    ap.add_argument(
        "-o",
        "--out",
        required=True,
        metavar="FILE",
        help="Output file path. If it ends with .gz, gzip output is produced.",
    )
    ap.add_argument(
        "inputs",
        nargs="+",
        metavar="IN",
        help="Input replicate files (2 or 3). Each line: chrom pos col3 col4. .gz supported.",
    )
    ap.epilog = (
        "Note: inputs must be sorted by (chrom,pos), otherwise the merge is not reliable.\n"
        "Example: bts-cx-merge -o merged.CX rep1.CX rep2.CX"
    )
    return ap


def merge_replicates(out_path: str, inputs: List[str]) -> None:
    if len(inputs) not in (2, 3):
        raise SystemExit("ERROR: please provide 2 or 3 input files.")

    hs = [open_text_read(p) for p in inputs]
    try:
        cur = [read_next(h) for h in hs]

        with open_text_write(out_path) as w:
            while True:
                kmin = min_key(cur)
                if kmin is None:
                    break

                out_chr: Optional[str] = None
                out_pos: Optional[int] = None
                s3 = 0
                s4 = 0

                for i in range(len(hs)):
                    if cur[i] is not None and cur[i][0] == kmin:
                        _, chrom, pos, v3, v4 = cur[i]
                        out_chr, out_pos = chrom, pos
                        s3 += v3
                        s4 += v4
                        cur[i] = read_next(hs[i])

                w.write(f"{out_chr}\t{out_pos}\t{s3}\t{s4}\n")

    finally:
        for h in hs:
            h.close()

    print(f"[{SUITE_NAME}] Done: {out_path}")


def main() -> None:
    args = build_parser().parse_args()
    merge_replicates(args.out, args.inputs)


if __name__ == "__main__":
    main()
