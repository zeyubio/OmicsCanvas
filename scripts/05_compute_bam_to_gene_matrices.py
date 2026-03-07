#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""OmicsCanvas: BAM -> gene-centric matrices (TSS / gene profile / TES)

Inputs
------
1) Sorted and indexed BAM file (.bam + .bai)
2) Gene BED file with 6 columns:
   chrom, start, end, feature_id, score, strand(+/-)

Outputs
-------
Three matrices with one row per feature_id:
  * <prefix>_tss_matrix.tsv
  * <prefix>_gene_profile_matrix.tsv
  * <prefix>_tes_matrix.tsv

Notes
-----
- This script normalizes coverage as:
    normalized = (A+C+G+T) / mapped_reads / read_length / lib_factor * norm_scale
  where lib_factor = 1 for single-end, 2 for paired-end.
- For negative-strand features, bins are reversed so all outputs are 5' -> 3'.
- Features too close to chromosome edges (window goes out of range) are skipped.

Author
------
OmicsCanvas project
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import re
import time
from typing import Any, Dict, Optional, Tuple, List

import numpy as np
import pandas as pd
import pysam


TYPES_READS: Dict[str, int] = {"single": 1, "paired": 2}


def load_config(path: str) -> Dict[str, Any]:
    """Load a JSON config file."""
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def cfg_get(cfg: Dict[str, Any], keys: List[str], default: Any) -> Any:
    """Nested dict getter."""
    cur: Any = cfg
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def validate_bins(distance: int,
                  tss_bins: int,
                  tes_bins: int,
                  gene_up_bins: int,
                  gene_body_bins: int,
                  gene_down_bins: int) -> None:
    """Validate divisibility constraints for binning."""
    if distance <= 0:
        raise ValueError("distance must be > 0")

    if tss_bins <= 0 or tes_bins <= 0:
        raise ValueError("tss_bins/tes_bins must be > 0")

    if gene_up_bins <= 0 or gene_body_bins <= 0 or gene_down_bins <= 0:
        raise ValueError("gene_*_bins must be > 0")

    if (2 * distance) % tss_bins != 0:
        raise ValueError(
            f"2*distance must be divisible by tss_bins. Got 2*{distance} % {tss_bins} != 0"
        )
    if (2 * distance) % tes_bins != 0:
        raise ValueError(
            f"2*distance must be divisible by tes_bins. Got 2*{distance} % {tes_bins} != 0"
        )
    if distance % gene_up_bins != 0:
        raise ValueError(
            f"distance must be divisible by gene_upstream_bins. Got {distance} % {gene_up_bins} != 0"
        )
    if distance % gene_down_bins != 0:
        raise ValueError(
            f"distance must be divisible by gene_downstream_bins. Got {distance} % {gene_down_bins} != 0"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="omicscanvas_bam_to_gene_matrices",
        description=(
            "Generate gene-centric coverage matrices around TSS/TES and across the gene body from a BAM file.\n"
            "Outputs three matrices: TSS, gene profile, TES."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--config",
        default=None,
        help=(
            "Optional JSON config file. Values in command line override config.\n"
            "Recommended for OmicsCanvas pipelines."
        ),
    )

    # Required inputs
    parser.add_argument(
        "-b",
        "--bam",
        required=True,
        help="Input BAM file (sorted and indexed, with .bai).",
    )
    parser.add_argument(
        "-g",
        "--bed",
        required=True,
        help=(
            "Input BED with 6 columns: chrom, start, end, feature_id, score, strand(+/-).\n"
            "Coordinates: 0-based, half-open [start, end)."
        ),
    )

    # Region definition
    parser.add_argument(
        "--distance",
        type=int,
        default=None,
        help="Flanking distance (bp) used for TSS/TES and gene upstream/downstream windows.",
    )

    # Bins (preferred names)
    parser.add_argument("--tss-bins", type=int, default=None, help="Number of bins in TSS ± distance window.")
    parser.add_argument("--tes-bins", type=int, default=None, help="Number of bins in TES ± distance window.")
    parser.add_argument(
        "--gene-up-bins",
        dest="gene_upstream_bins",
        type=int,
        default=None,
        help="Bins for upstream distance in gene profile.",
    )
    parser.add_argument("--gene-body-bins", type=int, default=None, help="Bins for scaled gene body in gene profile.")
    parser.add_argument(
        "--gene-down-bins",
        dest="gene_downstream_bins",
        type=int,
        default=None,
        help="Bins for downstream distance in gene profile.",
    )

    # Backward-compatible aliases (deprecated)
    parser.add_argument("--bins-start", dest="tss_bins", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--bins-end", dest="tes_bins", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--bins-gene-st", dest="gene_upstream_bins", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--bins-gene-body", dest="gene_body_bins", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--bins-gene-en", dest="gene_downstream_bins", type=int, default=None, help=argparse.SUPPRESS)

    # Legacy long-form names (keep parsing, hide from -h)
    parser.add_argument("--gene-upstream-bins", dest="gene_upstream_bins", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--gene-downstream-bins", dest="gene_downstream_bins", type=int, default=None, help=argparse.SUPPRESS)

    # Normalization
    parser.add_argument(
        "--read-length",
        type=int,
        default=150,
        help="Read length used for normalization.",
    )
    parser.add_argument(
        "--reads-type",
        choices=list(TYPES_READS.keys()),
        default="paired",
        help="Library type used for normalization: single or paired.",
    )
    parser.add_argument(
        "--norm-scale",
        type=float,
        default=1e8,
        help="Multiplicative scale factor applied after normalization.",
    )
    parser.add_argument(
        "--max-cov",
        type=float,
        default=500.0,
        help="Cap normalized coverage values at this threshold (winsorization).",
    )

    # Runtime
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=None,
        help="Number of worker processes.",
    )

    # Output
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. If set, output files will be written into this directory.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        help=(
            "Output prefix (sample ID). If outdir is provided, files will be written as outdir/<prefix>_*.\n"
            "If out-prefix contains a path separator, it will be used as-is."
        ),
    )
    parser.add_argument(
        "--naming",
        choices=["standard", "legacy"],
        default="standard",
        help=(
            "Output filename scheme.\n"
            "standard: *_tss_matrix.tsv, *_gene_profile_matrix.tsv, *_tes_matrix.tsv\n"
            "legacy: *_start_matrix.txt, *_gene_body_matrix.txt, *_end_matrix.txt"
        ),
    )
    parser.add_argument("--suffix-tss", default=None, help="Custom suffix for TSS matrix (overrides --naming).")
    parser.add_argument("--suffix-gene", default=None, help="Custom suffix for gene profile matrix (overrides --naming).")
    parser.add_argument("--suffix-tes", default=None, help="Custom suffix for TES matrix (overrides --naming).")

    return parser


def parse_args() -> argparse.Namespace:
    return build_parser().parse_args()


def _safe_prefix(out_prefix: str) -> str:
    """Strip common extensions if user accidentally includes them."""
    return re.sub(r"\.(txt|tsv|csv)$", "", out_prefix, flags=re.IGNORECASE)


def resolve_outputs(prefix: str,
                    naming: str,
                    suffix_tss: Optional[str],
                    suffix_gene: Optional[str],
                    suffix_tes: Optional[str]) -> Tuple[str, str, str]:
    if naming == "legacy":
        d_tss, d_gene, d_tes = "_start_matrix.txt", "_gene_body_matrix.txt", "_end_matrix.txt"
    else:
        d_tss, d_gene, d_tes = "_tss_matrix.tsv", "_gene_profile_matrix.tsv", "_tes_matrix.tsv"

    s_tss = suffix_tss if suffix_tss is not None else d_tss
    s_gene = suffix_gene if suffix_gene is not None else d_gene
    s_tes = suffix_tes if suffix_tes is not None else d_tes

    return f"{prefix}{s_gene}", f"{prefix}{s_tss}", f"{prefix}{s_tes}"


def run_start_bins(cov: np.ndarray,
                   distance: int,
                   tss_bins: int,
                   bed_sub: pd.DataFrame,
                   out_list: List[List[Any]]) -> None:
    """TSS matrix (TSS ± distance)."""
    step = int((2 * distance) / tss_bins)

    for st, en, fid, strand in bed_sub.values:
        if strand == "+":
            part = cov[st - distance: st + distance]
            if len(part) == 2 * distance:
                binned = part.reshape(-1, step).mean(axis=1)
                out_list.append([fid] + list(binned))
        else:
            part = cov[en - distance: en + distance]
            if len(part) == 2 * distance:
                binned = part.reshape(-1, step).mean(axis=1)
                out_list.append([fid] + list(binned)[::-1])


def run_end_bins(cov: np.ndarray,
                 distance: int,
                 tes_bins: int,
                 bed_sub: pd.DataFrame,
                 out_list: List[List[Any]]) -> None:
    """TES matrix (TES ± distance)."""
    step = int((2 * distance) / tes_bins)

    for st, en, fid, strand in bed_sub.values:
        if strand == "+":
            part = cov[en - distance: en + distance]
            if len(part) == 2 * distance:
                binned = part.reshape(-1, step).mean(axis=1)
                out_list.append([fid] + list(binned))
        else:
            part = cov[st - distance: st + distance]
            if len(part) == 2 * distance:
                binned = part.reshape(-1, step).mean(axis=1)
                out_list.append([fid] + list(binned)[::-1])


def run_gene_bins(cov: np.ndarray,
                  distance: int,
                  gene_up_bins: int,
                  gene_body_bins: int,
                  gene_down_bins: int,
                  bed_sub: pd.DataFrame,
                  out_list: List[List[Any]]) -> None:
    """Gene profile: upstream(distance) + scaled gene body + downstream(distance)."""
    up_step = int(distance / gene_up_bins)
    down_step = int(distance / gene_down_bins)

    for st, en, fid, strand in bed_sub.values:
        gene_len = int(en - st)
        if gene_len <= 0:
            continue

        if strand == "+":
            part_up = cov[st - distance: st]
            part_gene = cov[st: en]
            part_down = cov[en: en + distance]

            if len(part_up) != distance or len(part_down) != distance or len(part_gene) != gene_len:
                continue

            up = part_up.reshape(-1, up_step).mean(axis=1)
            down = part_down.reshape(-1, down_step).mean(axis=1)

            idx = np.floor(np.arange(gene_len) / float(gene_len) * gene_body_bins).astype(int)
            idx[idx >= gene_body_bins] = gene_body_bins - 1
            gene_df = pd.DataFrame({"idx": idx, "cov": part_gene})
            gene_mean = gene_df.groupby("idx")["cov"].mean().reindex(range(gene_body_bins), fill_value=np.nan)

            out_list.append([fid] + list(up) + list(gene_mean.values) + list(down))

        else:
            part_up = cov[en: en + distance]
            part_gene = cov[st: en]
            part_down = cov[st - distance: st]

            if len(part_up) != distance or len(part_down) != distance or len(part_gene) != gene_len:
                continue

            # for negative strand, reverse to keep 5'->3'
            up = part_up.reshape(-1, down_step).mean(axis=1)[::-1]
            down = part_down.reshape(-1, up_step).mean(axis=1)[::-1]

            idx = np.floor(np.arange(gene_len) / float(gene_len) * gene_body_bins).astype(int)
            idx[idx >= gene_body_bins] = gene_body_bins - 1
            gene_df = pd.DataFrame({"idx": idx, "cov": part_gene})
            gene_mean = gene_df.groupby("idx")["cov"].mean().reindex(range(gene_body_bins), fill_value=np.nan)

            out_list.append([fid] + list(up) + list(gene_mean.values[::-1]) + list(down))


def run_one_chrom(chrom: str, args_dict: Dict[str, Any], gene_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Worker: compute matrices for one chromosome."""
    distance = int(args_dict["distance"])
    tss_bins = int(args_dict["tss_bins"])
    tes_bins = int(args_dict["tes_bins"])
    gene_up_bins = int(args_dict["gene_up_bins"])
    gene_body_bins = int(args_dict["gene_body_bins"])
    gene_down_bins = int(args_dict["gene_down_bins"])

    local_time = time.localtime()
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', local_time)} : processing {chrom}")

    ch_gene = gene_df.loc[chrom]
    if isinstance(ch_gene, pd.Series):
        ch_gene = ch_gene.to_frame().T

    with pysam.AlignmentFile(args_dict["bam"], "rb") as bam:
        mapped = bam.mapped
        if mapped == 0:
            raise RuntimeError(f"No mapped reads found in BAM: {args_dict['bam']}")
        ch_len = bam.get_reference_length(chrom)
        map_A, map_C, map_G, map_T = bam.count_coverage(contig=chrom, start=0, stop=ch_len)

    cov_raw = np.array(map_A) + np.array(map_C) + np.array(map_G) + np.array(map_T)
    cov_norm = (cov_raw / float(mapped) / float(args_dict["read_length"]) / float(TYPES_READS[args_dict["reads_type"]]) * float(args_dict["norm_scale"]))
    cov_norm = cov_norm.astype(float)
    cov_norm[cov_norm > float(args_dict["max_cov"])] = float(args_dict["max_cov"])

    out_tss: List[List[Any]] = []
    out_tes: List[List[Any]] = []
    out_gene: List[List[Any]] = []

    bed_sub = ch_gene[["st", "en", "ID", "zf"]]

    run_start_bins(cov_norm, distance, tss_bins, bed_sub, out_tss)
    run_end_bins(cov_norm, distance, tes_bins, bed_sub, out_tes)
    run_gene_bins(cov_norm, distance, gene_up_bins, gene_body_bins, gene_down_bins, bed_sub, out_gene)

    gene_cols = ["ID"] + [f"bin_{i+1}" for i in range(gene_up_bins + gene_body_bins + gene_down_bins)]
    tss_cols = ["ID"] + [f"bin_{i+1}" for i in range(tss_bins)]
    tes_cols = ["ID"] + [f"bin_{i+1}" for i in range(tes_bins)]

    df_gene = pd.DataFrame(out_gene, columns=gene_cols)
    df_tss = pd.DataFrame(out_tss, columns=tss_cols)
    df_tes = pd.DataFrame(out_tes, columns=tes_cols)

    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', local_time)} : finished {chrom}")
    return df_gene, df_tss, df_tes


def main() -> None:
    args = parse_args()

    cfg: Dict[str, Any] = {}
    if args.config:
        cfg = load_config(args.config)

    # Resolve parameters (CLI overrides config)
    distance = args.distance if args.distance is not None else int(cfg_get(cfg, ["genomic_bins", "distance"], 2000))

    tss_bins = args.tss_bins if args.tss_bins is not None else int(cfg_get(cfg, ["genomic_bins", "tss_bins"], 100))
    tes_bins = args.tes_bins if args.tes_bins is not None else int(cfg_get(cfg, ["genomic_bins", "tes_bins"], 100))

    gene_up_bins = args.gene_upstream_bins if args.gene_upstream_bins is not None else int(cfg_get(cfg, ["genomic_bins", "gene_body", "upstream_bins"], 50))
    gene_body_bins = args.gene_body_bins if args.gene_body_bins is not None else int(cfg_get(cfg, ["genomic_bins", "gene_body", "body_bins"], 100))
    gene_down_bins = args.gene_downstream_bins if args.gene_downstream_bins is not None else int(cfg_get(cfg, ["genomic_bins", "gene_body", "downstream_bins"], 50))

    threads = args.threads if args.threads is not None else int(cfg_get(cfg, ["threads"], 5))

    outdir_cfg = cfg_get(cfg, ["paths", "matrix_dir"], None)
    outdir = args.outdir if args.outdir is not None else outdir_cfg

    validate_bins(distance, tss_bins, tes_bins, gene_up_bins, gene_body_bins, gene_down_bins)

    # Build prefix path
    prefix = _safe_prefix(args.out_prefix)
    if os.path.sep in prefix:
        prefix_path = prefix
    else:
        if outdir:
            os.makedirs(outdir, exist_ok=True)
            prefix_path = os.path.join(outdir, prefix)
        else:
            prefix_path = prefix

    out_gene, out_tss, out_tes = resolve_outputs(
        prefix_path,
        args.naming,
        args.suffix_tss,
        args.suffix_gene,
        args.suffix_tes,
    )

    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    print(f"BAM: {args.bam}")
    print(f"BED: {args.bed}")
    print(f"distance={distance}; tss_bins={tss_bins}; gene_bins={gene_up_bins}+{gene_body_bins}+{gene_down_bins}; tes_bins={tes_bins}")
    print(f"threads={threads}")
    print(f"Outputs:\n  {out_tss}\n  {out_gene}\n  {out_tes}")

    bed = pd.read_csv(
        args.bed,
        sep="\t",
        header=None,
        names=["ch", "st", "en", "ID", "score", "zf"],
        usecols=[0, 1, 2, 3, 4, 5],
    )
    gene = bed[["ch", "st", "en", "ID", "zf"]].copy()
    gene["st"] = pd.to_numeric(gene["st"], errors="coerce").astype(int)
    gene["en"] = pd.to_numeric(gene["en"], errors="coerce").astype(int)
    gene = gene.dropna().set_index("ch")

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        chroms_bam = list(bam.references)

    chroms_bed = sorted(list(gene.index.unique()))
    chroms = list(np.intersect1d(chroms_bed, chroms_bam))

    if len(chroms) == 0:
        raise RuntimeError("No common chromosomes found between BED and BAM.")

    args_dict: Dict[str, Any] = {
        "bam": args.bam,
        "distance": distance,
        "tss_bins": tss_bins,
        "tes_bins": tes_bins,
        "gene_up_bins": gene_up_bins,
        "gene_body_bins": gene_body_bins,
        "gene_down_bins": gene_down_bins,
        "read_length": args.read_length,
        "reads_type": args.reads_type,
        "norm_scale": args.norm_scale,
        "max_cov": args.max_cov,
    }

    with mp.Pool(processes=threads) as pool:
        results = pool.starmap(run_one_chrom, [(ch, args_dict, gene) for ch in chroms])

    df_gene_list = [r[0] for r in results]
    df_tss_list = [r[1] for r in results]
    df_tes_list = [r[2] for r in results]

    all_gene = pd.concat(df_gene_list, ignore_index=True)
    all_tss = pd.concat(df_tss_list, ignore_index=True)
    all_tes = pd.concat(df_tes_list, ignore_index=True)

    all_tss.to_csv(out_tss, sep="\t", index=False)
    all_gene.to_csv(out_gene, sep="\t", index=False)
    all_tes.to_csv(out_tes, sep="\t", index=False)

    print(f"End time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
    print("Done.")


if __name__ == "__main__":
    main()
