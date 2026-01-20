#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CXGeneMatrix

Build feature-centric methylation matrices from a single-context CX file (CG/CHG/CHH).
The default input is a gene BED (6 columns), but any interval set works (e.g., TE BED),
as long as the BED follows the required format.

Method
- Extend each feature by "distance" bp on both sides.
- Discretize upstream / feature body / downstream into bins.
- Scan CX sites (pos, me, al) and assign each site to the corresponding bin.
- Outputs:
  1) per-feature bin tables (upstream / body / downstream)
  2) a global methylation profile (ratio = me/al per bin, concatenated upstream→body→downstream)

Input formats
1) BED (recommended 6 columns; strand is required)
   chrom  start  end  id  score  strand
   - Use the same coordinate convention as your CX file (BED is often 0-based while CX may be 1-based).
     If your pipeline already uses consistent coordinates, no conversion is needed.
2) CX (4 columns, no header)
   ch  pos  me  al
   - ch: chromosome
   - pos: genomic position (integer)
   - me: methylated read/support count
   - al: total coverage/support count

Outputs
- <out_dir>/<sample>_<context>_upstream_bins<BINS_UP>.tsv
- <out_dir>/<sample>_<context>_body_bins<BINS_BODY>.tsv
- <out_dir>/<sample>_<context>_downstream_bins<BINS_DOWN>.tsv
- <out_dir>/<sample>_<context>_profile.tsv
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Dict, List

import numpy as np
import pandas as pd


# ----------------------------
# CLI
# ----------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="CXGeneMatrix",
        description="Generate feature-centric methylation matrices and a global profile from a single-context CX file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s", "--sample", required=True,
        help="Sample/tissue name used as output prefix (e.g., SRR9321764 or S2).",
    )
    parser.add_argument(
        "-c", "--context", default="CHH", choices=["CG", "CHG", "CHH"],
        help="Methylation context (CG/CHG/CHH).",
    )

    # BED input
    parser.add_argument(
        "-b", "--bed", required=True,
        help=(
            "BED path (at least 6 columns): chrom, start, end, id, score, strand.\n"
            "The strand column must be '+' or '-'."
        ),
    )

    # CX input
    parser.add_argument(
        "--cx-file", default=None,
        help=(
            "CX path (4 columns): ch, pos, me, al.\n"
            "If omitted, the default is: <cx-dir>/<sample>_<context>.CX"
        ),
    )
    parser.add_argument(
        "--cx-dir", default="meth_data",
        help="Directory holding CX files (used to infer the default CX path).",
    )

    # Binning
    parser.add_argument(
        "--distance", type=int, default=2000,
        help="Upstream/downstream extension length in bp.",
    )
    parser.add_argument(
        "--bins-up", type=int, default=50,
        help="Number of bins for the upstream region.",
    )
    parser.add_argument(
        "--bins-body", type=int, default=100,
        help="Number of bins for the feature body.",
    )
    parser.add_argument(
        "--bins-down", type=int, default=50,
        help="Number of bins for the downstream region.",
    )

    # Chrom filter
    parser.add_argument(
        "--chrom-prefix", default="Chr",
        help="Keep only records whose chromosome starts with this prefix (set to empty to disable).",
    )

    # Output
    parser.add_argument(
        "-o", "--out-dir", default="CX_gene",
        help="Output directory.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Overwrite existing outputs (default: fail if outputs already exist).",
    )

    parser.epilog = (
        "Example:\n"
        "  python cx_gene_matrix.en.py \\\n"
        "    -s SRR9321764 -c CHH \\\n"
        "    -b genes.bed \\\n"
        "    --cx-dir meth_data \\\n"
        "    --distance 2000 --bins-up 50 --bins-body 100 --bins-down 50 \\\n"
        "    -o CX_gene\n"
    )
    return parser


def parse_args(argv: List[str]) -> argparse.Namespace:
    return build_parser().parse_args(argv)


# ----------------------------
# Core logic (kept close to the original)
# ----------------------------

def _select_non_overlapping(intervals: np.ndarray,
                            kept: List[List],
                            rejected: List[List]) -> None:
    """Greedy selection of non-overlapping intervals (intervals must be sorted by start)."""
    kept.append(list(intervals[0]))
    for i in range(1, len(intervals)):
        if intervals[i][0] >= kept[-1][1]:
            kept.append(list(intervals[i]))
        else:
            rejected.append(list(intervals[i]))


def split_overlapping_features(intervals: np.ndarray, max_rounds: int = 10) -> Dict[int, List[List]]:
    """
    Split potentially overlapping features into multiple groups where features are non-overlapping within each group.
    This allows a linear pointer scan of CX sites for each group.
    """
    groups: Dict[int, List[List]] = {"trash": []}  # type: ignore[assignment]

    for r in range(max_rounds):
        groups[r] = []
        if r == 0:
            _select_non_overlapping(intervals, groups[r], groups["trash"])  # type: ignore[index]
        else:
            if len(groups["trash"]) == 0:  # type: ignore[arg-type]
                del groups[r]
                del groups["trash"]  # type: ignore[arg-type]
                break
            tmp = groups["trash"]  # type: ignore[assignment]
            groups["trash"] = []   # type: ignore[assignment]
            _select_non_overlapping(np.asarray(tmp, dtype=object), groups[r], groups["trash"])  # type: ignore[index]

    if "trash" in groups:
        del groups["trash"]  # type: ignore[arg-type]
    return groups  # type: ignore[return-value]


def scan_one_chrom(features: List[List],
                   sites: np.ndarray,
                   distance: int,
                   bins_up: int,
                   bins_body: int,
                   bins_down: int,
                   out_up: List[List],
                   out_body: List[List],
                   out_down: List[List]) -> None:
    """
    Scan one chromosome and assign CX sites to upstream/body/downstream bins for each feature.

    sites: shape (N, 3) -> [pos, me, al]
    features: [ [po1, po2, strand, id], ... ] (po1/po2 already include flanks)
    """
    pointer = 0

    for (po1, po2, strand, fid) in features:
        po1 = int(po1)
        po2 = int(po2)

        if strand == "+":
            for j in range(pointer, len(sites)):
                pos = int(sites[j][0])
                me = sites[j][1]
                al = sites[j][2]

                if po1 < pos < po1 + distance:
                    b = int(np.ceil((pos - po1 - 1) / distance * bins_up))
                    out_up.append([fid, b, me, al])

                if po1 + distance <= pos <= po2 - distance:
                    body_len = (po2 - po1 - 2 * distance + 1)
                    step = np.ceil(body_len / bins_body) if body_len > 0 else 1
                    b = int(np.ceil((pos - po1 - distance) / step))
                    out_body.append([fid, b, me, al])

                if po2 - distance < pos < po2:
                    b = int(np.ceil((pos - po2 + distance - 1) / (distance / bins_down)))
                    out_down.append([fid, b, me, al])

                if pos >= po2:
                    pointer = j
                    break

        else:  # strand == "-"
            for j in range(pointer, len(sites)):
                pos = int(sites[j][0])
                me = sites[j][1]
                al = sites[j][2]

                if po1 < pos < po1 + distance:
                    b = int((po1 + distance - pos - 1) / (distance / bins_down))
                    out_down.append([fid, b, me, al])

                if po1 + distance <= pos <= po2 - distance:
                    body_len = (po2 - po1 - 2 * distance + 1)
                    step = np.ceil(body_len / bins_body) if body_len > 0 else 1
                    b = int((po2 - distance - pos) / step)
                    out_body.append([fid, b, me, al])

                if po2 - distance < pos < po2:
                    b = int(np.ceil((po2 - pos - 1) / (distance / bins_up)))
                    out_up.append([fid, b, me, al])

                if pos >= po2:
                    pointer = j
                    break


# ----------------------------
# I/O helpers and pipeline
# ----------------------------

def _require_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{label} not found: {path}")


def _safe_mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _check_overwrite(paths: List[str], overwrite: bool) -> None:
    existing = [p for p in paths if os.path.exists(p)]
    if existing and not overwrite:
        msg = "Output files already exist (use --overwrite to replace):\n  " + "\n  ".join(existing)
        raise FileExistsError(msg)
    if overwrite:
        for p in existing:
            os.remove(p)


def load_bed(bed_path: str, chrom_prefix: str, distance: int) -> pd.DataFrame:
    bed = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        comment="#",
        dtype={0: str},
    )

    if bed.shape[1] < 6:
        raise ValueError("BED must have at least 6 columns: chrom, start, end, id, score, strand.")

    bed = bed.iloc[:, :6].copy()
    bed.columns = ["ch", "st", "en", "id", "score", "strand"]

    if chrom_prefix:
        bed = bed[bed["ch"].astype(str).str.startswith(chrom_prefix)].copy()

    if not set(bed["strand"].unique()).issubset({"+", "-"}):
        raise ValueError("The 6th BED column (strand) must be '+' or '-'.")

    bed["st"] = bed["st"].astype(int)
    bed["en"] = bed["en"].astype(int)
    bed["po1"] = bed["st"] - distance
    bed["po2"] = bed["en"] + distance
    bed.sort_values(["ch", "st"], inplace=True)
    return bed


def load_cx(cx_path: str, chrom_prefix: str) -> pd.DataFrame:
    cx = pd.read_csv(
        cx_path,
        sep="\t",
        header=None,
        names=["ch", "pos", "me", "al"],
        dtype={"ch": str, "pos": int, "me": float, "al": float},
    )
    if chrom_prefix:
        cx = cx[cx["ch"].astype(str).str.startswith(chrom_prefix)].copy()
    cx.sort_values(["ch", "pos"], inplace=True)
    return cx


def write_profile(up_df: pd.DataFrame,
                  body_df: pd.DataFrame,
                  down_df: pd.DataFrame,
                  bins_up: int,
                  bins_body: int,
                  bins_down: int,
                  out_path: str) -> None:
    def _agg(df: pd.DataFrame, n_bins: int) -> pd.DataFrame:
        agg = df.groupby("bin")[["me", "al"]].sum()
        agg = agg.reindex(range(1, n_bins + 1), fill_value=0)
        ratio = np.where(agg["al"].to_numpy() > 0, agg["me"].to_numpy() / agg["al"].to_numpy(), np.nan)
        agg["ratio"] = ratio
        return agg

    a = _agg(up_df, bins_up)["ratio"].tolist()
    b = _agg(body_df, bins_body)["ratio"].tolist()
    c = _agg(down_df, bins_down)["ratio"].tolist()

    profile = pd.DataFrame({"ratio": a + b + c})
    profile.to_csv(out_path, sep="\t", index=False)


def main(argv: List[str]) -> int:
    args = parse_args(argv)

    cx_path = args.cx_file or os.path.join(args.cx_dir, f"{args.sample}_{args.context}.CX")

    _require_file(args.bed, "BED")
    _require_file(cx_path, "CX")
    _safe_mkdir(args.out_dir)

    prefix = f"{args.sample}_{args.context}"
    out_up = os.path.join(args.out_dir, f"{prefix}_upstream_bins{args.bins_up}.tsv")
    out_body = os.path.join(args.out_dir, f"{prefix}_body_bins{args.bins_body}.tsv")
    out_down = os.path.join(args.out_dir, f"{prefix}_downstream_bins{args.bins_down}.tsv")
    out_profile = os.path.join(args.out_dir, f"{prefix}_profile.tsv")

    _check_overwrite([out_up, out_body, out_down, out_profile], args.overwrite)

    sys.stderr.write(f"[INFO] start: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\n")
    sys.stderr.write(f"[INFO] BED: {args.bed}\n")
    sys.stderr.write(f"[INFO] CX : {cx_path}\n")
    sys.stderr.write(
        f"[INFO] distance={args.distance}, bins_up={args.bins_up}, bins_body={args.bins_body}, bins_down={args.bins_down}\n"
    )

    bed = load_bed(args.bed, args.chrom_prefix, args.distance)
    cx = load_cx(cx_path, args.chrom_prefix)

    chroms = sorted(bed["ch"].unique().tolist())
    if len(chroms) == 0:
        raise ValueError("No chromosomes left after filtering. Check --chrom-prefix.")

    up_rows: List[List] = []
    body_rows: List[List] = []
    down_rows: List[List] = []

    for ch in chroms:
        feats = bed[bed["ch"] == ch][["po1", "po2", "strand", "id"]].to_numpy(dtype=object)
        if feats.shape[0] == 0:
            continue

        sites = cx[cx["ch"] == ch][["pos", "me", "al"]].to_numpy()
        if sites.shape[0] == 0:
            sys.stderr.write(f"[WARN] {ch}: no CX records; skipped\n")
            continue

        groups = split_overlapping_features(feats)
        for k in sorted(groups.keys()):
            scan_one_chrom(
                groups[k],
                sites,
                args.distance,
                args.bins_up,
                args.bins_body,
                args.bins_down,
                up_rows,
                body_rows,
                down_rows,
            )
        sys.stderr.write(f"[INFO] {ch} done: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\n")

    up_df = pd.DataFrame(up_rows, columns=["id", "bin", "me", "al"])
    body_df = pd.DataFrame(body_rows, columns=["id", "bin", "me", "al"])
    down_df = pd.DataFrame(down_rows, columns=["id", "bin", "me", "al"])

    up_df.to_csv(out_up, sep="\t", index=False)
    body_df.to_csv(out_body, sep="\t", index=False)
    down_df.to_csv(out_down, sep="\t", index=False)

    write_profile(up_df, body_df, down_df, args.bins_up, args.bins_body, args.bins_down, out_profile)

    sys.stderr.write(f"[INFO] done: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\n")
    sys.stderr.write(f"[INFO] outputs:\n  {out_up}\n  {out_body}\n  {out_down}\n  {out_profile}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
