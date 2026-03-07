#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import pandas as pd
pysam = None  # imported lazily in main()



# ----------------------------
# Filters
# ----------------------------
def pass_filters(r, min_mapq=0, include_duplicates=False, include_qcfail=False):
    """
    Return True if read passes strict filters.
    Strict by default:
      - mapped
      - not secondary
      - not supplementary
      - (optional) not duplicate
      - (optional) not QC fail
      - MAPQ >= min_mapq
    """
    if r.is_unmapped:
        return False
    if r.is_secondary or r.is_supplementary:
        return False
    if (not include_duplicates) and r.is_duplicate:
        return False
    if (not include_qcfail) and r.is_qcfail:
        return False
    if r.mapping_quality < min_mapq:
        return False
    return True


def is_valid_proper_pair(r):
    """
    Require proper paired and both mates mapped.
    For strict PE fragment counting we default to proper pairs only.
    """
    if not r.is_paired:
        return False
    if r.mate_is_unmapped:
        return False
    if r.is_unmapped:
        return False
    if not r.is_proper_pair:
        return False
    return True


# ----------------------------
# IO
# ----------------------------
def read_bed(path, one_based=False):
    """
    BED expected columns (at least): chr, start, end, ID
    pysam uses 0-based, end-exclusive coordinates.

    If your "BED" was actually 1-based inclusive (common mistake), use --bed-one-based
    """
    bed = pd.read_csv(
        path, sep="\t", header=None, comment="#",
        usecols=[0, 1, 2, 3],
        names=["ch", "st", "en", "ID"],
        dtype={"ch": str, "st": int, "en": int, "ID": str},
    ).dropna()

    if one_based:
        # Convert 1-based inclusive to 0-based half-open:
        # [st, en] -> [st-1, en)
        bed["st"] = bed["st"] - 1
        # en stays the same as half-open end

    return bed


def read_length(path):
    """
    length table: 2 columns, no header: ID <tab> length_bp
    """
    df = pd.read_csv(
        path, sep="\t", header=None, comment="#",
        names=["ID", "length_bp"],
        dtype={"ID": str, "length_bp": float},
    ).dropna()

    df = df[df["length_bp"] > 0].copy()
    return df.set_index("ID")["length_bp"]


# ----------------------------
# Counting
# ----------------------------
def count_libsize_se(bam, min_mapq, include_duplicates, include_qcfail):
    """
    Strict library size for SE: count reads across whole BAM using same filters.
    """
    n = 0
    for r in bam.fetch(until_eof=True):
        if pass_filters(r, min_mapq=min_mapq,
                        include_duplicates=include_duplicates,
                        include_qcfail=include_qcfail):
            n += 1
    return float(n)


def count_libsize_pe(bam, min_mapq, include_duplicates, include_qcfail, require_proper=True):
    """
    Strict library size for PE: count fragments across whole BAM.

    Strategy:
      - count only read1 to count each fragment once
      - require proper pair by default (recommended strict mode)
    """
    n = 0
    for r in bam.fetch(until_eof=True):
        if not r.is_paired:
            continue
        if require_proper and (not is_valid_proper_pair(r)):
            continue
        if not r.is_read1:
            continue
        if pass_filters(r, min_mapq=min_mapq,
                        include_duplicates=include_duplicates,
                        include_qcfail=include_qcfail):
            n += 1
    return float(n)


def count_interval_se(bam, ch, st, en, min_mapq, include_duplicates, include_qcfail):
    """
    SE interval count: reads overlapping interval.
    """
    c = 0
    for r in bam.fetch(ch, st, en):
        if pass_filters(r, min_mapq=min_mapq,
                        include_duplicates=include_duplicates,
                        include_qcfail=include_qcfail):
            c += 1
    return float(c)


def count_interval_pe_fragments(bam, ch, st, en, min_mapq, include_duplicates, include_qcfail, require_proper=True):
    """
    PE strict interval count: count FRAGMENTS (templates) overlapping interval.
    Deduplicate by query_name so a fragment is counted once even if both mates overlap.

    Rule (strict):
      - require proper pair by default
      - include a fragment if either mate overlaps interval (because fetch sees whichever overlaps)
      - exclude secondary/supplementary/unmapped/duplicate/qcfail per filters
    """
    qnames = set()
    for r in bam.fetch(ch, st, en):
        if not r.is_paired:
            continue
        if require_proper and (not is_valid_proper_pair(r)):
            continue
        if pass_filters(r, min_mapq=min_mapq,
                        include_duplicates=include_duplicates,
                        include_qcfail=include_qcfail):
            qnames.add(r.query_name)
    return float(len(qnames))


# ----------------------------
# Main
# ----------------------------
def build_parser():
    """Build CLI parser (English help strings)."""
    p = argparse.ArgumentParser(
        prog="omicscanvas_bam_to_fpkm.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Compute strict FPKM/TPM from a BAM using pysam (single-end or paired-end).

Paired-end mode counts FRAGMENTS (deduplicated by query_name within each interval) and uses a consistent library-size definition.""",
        epilog="""Examples:
  # single-end
  python omicscanvas_bam_to_fpkm.py -b genes.bed -l gene_len.txt -m sample.bam -o out.tsv --mode se

  # paired-end (strict fragments)
  python omicscanvas_bam_to_fpkm.py -b genes.bed -l gene_len.txt -m sample.bam -o out.tsv --mode pe --min-mapq 10

Notes:
  - BAM must be coordinate-sorted and indexed (.bai).
  - BED must be 0-based, end-exclusive. If your BED is 1-based inclusive, use --bed-one-based.
""",
    )

    p.add_argument("-b", "--bed", required=True,
                   help="BED with at least 4 columns: chrom, start, end, gene_id (0-based, end-exclusive).")
    p.add_argument("-l", "--length", required=True,
                   help="Gene length table (2 columns, no header): gene_id<TAB>length_bp.")
    p.add_argument("-m", "--bam", required=True,
                   help="Input BAM (coordinate-sorted and indexed).")
    p.add_argument("-o", "--out", required=True,
                   help="Output TSV (index=gene_id; columns: FPKM, TPM, counts, length_bp).")

    p.add_argument("--mode", choices=["auto", "se", "pe"], default="auto",
                   help="Read layout: auto-detect from BAM, or force se/pe.")

    p.add_argument("--min-mapq", type=int, default=0,
                   help="Minimum MAPQ (mapping quality) to keep a read.")
    p.add_argument("--include-duplicates", action="store_true",
                   help="Include PCR/optical duplicates (default: excluded).")
    p.add_argument("--include-qcfail", action="store_true",
                   help="Include QC-failed reads (default: excluded).")

    p.add_argument("--require-proper-pair", action="store_true", default=True,
                   help="(PE) Require proper pairs (recommended).")
    p.add_argument("--allow-improper-pair", action="store_true",
                   help="(PE) Allow improper pairs (disables proper-pair requirement).")

    p.add_argument("--bed-one-based", action="store_true",
                   help="Convert a 1-based inclusive BED to 0-based half-open coordinates.")

    p.add_argument("--threads", type=int, default=1,
                   help="Threads for BAM BGZF decompression in pysam.")
    p.add_argument("--progress", type=int, default=1000,
                   help="Print progress every N intervals (0 disables).")

    return p


def detect_mode(bam):
    """
    Detect whether BAM is paired-end by inspecting first few reads.
    """
    n = 0
    paired = False
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped:
            continue
        paired = r.is_paired
        n += 1
        if n >= 200:
            break
    return "pe" if paired else "se"


def main():
    parser = build_parser()
    args = parser.parse_args()

    # Lazy import so `-h/--help` works even if pysam is not installed
    global pysam
    try:
        import pysam as _pysam
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError("pysam is required. Install: pip install pysam") from e
    pysam = _pysam


    require_proper = args.require_proper_pair and (not args.allow_improper_pair)

    # Open BAM with threads
    bam = pysam.AlignmentFile(args.bam, "rb", threads=max(1, args.threads))
    try:
        bam.check_index()
    except Exception:
        raise RuntimeError("BAM index not found or invalid. Please run: samtools index <bam>")

    mode = args.mode
    if mode == "auto":
        mode = detect_mode(bam)

    # After detect_mode, reopen BAM because until_eof consumed iterator
    bam.close()
    bam = pysam.AlignmentFile(args.bam, "rb", threads=max(1, args.threads))

    bed = read_bed(args.bed, one_based=args.bed_one_based)
    length_bp = read_length(args.length)

    # Align lengths to BED IDs
    # (If missing lengths exist, drop those IDs)
    bed_ids = bed["ID"].astype(str)
    missing = (~bed_ids.isin(length_bp.index)).sum()
    if missing > 0:
        print(f"[WARN] {missing} BED IDs missing in length table; they will be dropped.", file=sys.stderr)
        bed = bed[bed["ID"].isin(length_bp.index)].copy()

    # Library size N (strict, same filters)
    print(f"[INFO] Mode: {mode}", file=sys.stderr)
    print("[INFO] Counting library size (strict)...", file=sys.stderr)

    # Important: until_eof reads entire BAM, then we reopen for interval counting
    if mode == "se":
        N = count_libsize_se(
            bam, min_mapq=args.min_mapq,
            include_duplicates=args.include_duplicates,
            include_qcfail=args.include_qcfail
        )
    else:
        N = count_libsize_pe(
            bam, min_mapq=args.min_mapq,
            include_duplicates=args.include_duplicates,
            include_qcfail=args.include_qcfail,
            require_proper=require_proper
        )

    if N <= 0:
        raise RuntimeError("Library size N <= 0 after filtering. Check BAM/filter settings.")

    print(f"[INFO] Library size N = {N}", file=sys.stderr)

    # Reopen BAM for interval counting (because we consumed it during libsize scan)
    bam.close()
    bam = pysam.AlignmentFile(args.bam, "rb", threads=max(1, args.threads))

    # Count per interval
    counts = []
    for i, (ch, st, en, gid) in enumerate(bed[["ch", "st", "en", "ID"]].itertuples(index=False, name=None), start=1):
        if args.progress and (i % args.progress == 0):
            print(f"[INFO] Processed {i} intervals...", file=sys.stderr)

        if mode == "se":
            c = count_interval_se(
                bam, ch, int(st), int(en),
                min_mapq=args.min_mapq,
                include_duplicates=args.include_duplicates,
                include_qcfail=args.include_qcfail
            )
        else:
            c = count_interval_pe_fragments(
                bam, ch, int(st), int(en),
                min_mapq=args.min_mapq,
                include_duplicates=args.include_duplicates,
                include_qcfail=args.include_qcfail,
                require_proper=require_proper
            )
        counts.append((gid, c))

    df = pd.DataFrame(counts, columns=["ID", "counts"]).set_index("ID")
    df["length_bp"] = length_bp

    # Drop any remaining missing length
    df = df.dropna(subset=["length_bp"]).copy()

    # FPKM = 1e9 * counts / (length_bp * N)
    df["FPKM"] = (1e9 * df["counts"]) / (df["length_bp"] * N)

    # TPM
    df["RPK"] = df["counts"] / (df["length_bp"] / 1000.0)
    sum_rpk = df["RPK"].sum()
    if sum_rpk <= 0:
        raise RuntimeError("Sum of RPK <= 0, cannot compute TPM. Check counts/length.")
    df["TPM"] = 1e6 * df["RPK"] / sum_rpk


    # Output
    df = df[["FPKM", "TPM", "counts", "length_bp"]].sort_index()
    df.to_csv(args.out, sep="\t", index=True)

    print(f"[DONE] Output written: {args.out}", file=sys.stderr)
    # print(
    #     "[NOTE] PE mode counts fragments by query_name dedup within each interval. "
    #     "If genes overlap, the same fragment may be counted for multiple genes.\n"
    #     "双端模式按fragment（query_name）去重；若基因区间重叠，同一fragment可能被多个基因计入。",
    #     file=sys.stderr
    # )


if __name__ == "__main__":
    main()
