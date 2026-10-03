#!/usr/bin/env python3
"""Genotype the Chr23997 intron-2 PAV directly from BAM split-read evidence.

The 18-species assemblies show several overlapping, but non-identical,
deletions in intron 2.  This script therefore scores a *core-loss state* rather
than forcing every sample onto the 221-bp assembly allele.  A DEL-supporting
template must contain two high-quality Chr4 alignment segments whose reference
gap spans the pre-defined intron core.  Continuous alignments spanning the
same core provide reference-allele evidence.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import pysam


CORE_START = 88980879
CORE_END = 88981002


def gap_spans_core(left_end: int, right_start: int, core_start: int, core_end: int,
                   tolerance: int = 20, min_gap: int = 80) -> bool:
    """Whether an ordered reference gap is a candidate deletion of the core."""
    return (
        right_start - left_end >= min_gap
        and left_end <= core_start + tolerance
        and right_start >= core_end - tolerance
    )


def continuous_core_span(blocks, core_start: int, core_end: int, anchor: int = 20) -> bool:
    """Whether one aligned reference block continuously spans the core."""
    return any(start <= core_start - anchor and end >= core_end + anchor for start, end in blocks)


def call_genotype(flank_depth: float, core_depth: float, ref_templates: int,
                  del_templates: int) -> str:
    """Conservative genotype based on breakpoint and continuous-span evidence."""
    if flank_depth < 8:
        return "./."
    ratio = core_depth / flank_depth
    # A core-spanning CIGAR block is useful but too stringent for interspecific
    # short reads: ordinary small indels split many otherwise intact alignments.
    # Continuous depth through the core is therefore accepted as reference
    # evidence only when no deletion-breakpoint template is observed.
    if del_templates >= 3 and ref_templates < 2 and ratio <= 0.30:
        return "1/1"
    if del_templates >= 3 and (ref_templates >= 3 or ratio >= 0.50):
        return "0/1"
    if del_templates == 0 and (ref_templates >= 3 or ratio >= 0.50):
        return "0/0"
    return "./."


def _cigar_blocks(read):
    """Return aligned reference blocks, splitting at deletions and skipped bases."""
    ref_pos = read.reference_start
    block_start = ref_pos
    blocks = []
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):  # M, =, X
            ref_pos += length
        elif op in (2, 3):  # D, N
            if block_start < ref_pos:
                blocks.append((block_start, ref_pos))
            ref_pos += length
            block_start = ref_pos
        elif op in (1, 4, 5, 6):  # I, soft/hard clip, pad
            continue
    if block_start < ref_pos:
        blocks.append((block_start, ref_pos))
    return blocks


def _sa_segments(read, chrom: str, min_mapq: int):
    """Yield primary/supplementary and SA-tag segments on the target chromosome."""
    segments = [(read.reference_start, read.reference_end, read.mapping_quality)]
    if not read.has_tag("SA"):
        return segments
    for raw in read.get_tag("SA").strip(";").split(";"):
        fields = raw.split(",")
        if len(fields) < 6 or fields[0] != chrom or int(fields[4]) < min_mapq:
            continue
        start = int(fields[1]) - 1
        ref_length = sum(
            int(length) for length, op in re.findall(r"(\d+)([MDN=X])", fields[3])
        )
        segments.append((start, start + ref_length, int(fields[4])))
    return segments


def mean_depth(bam, chrom: str, start: int, end: int) -> float:
    depths = []
    for column in bam.pileup(chrom, start, end, truncate=True, stepper="samtools", min_mapping_quality=20):
        if start <= column.reference_pos < end:
            depths.append(column.nsegments)
    return sum(depths) / len(depths) if depths else 0.0


def collect_evidence(bam_path: str, chrom: str, core_start: int, core_end: int,
                     scan_flank: int = 1500):
    """Collect unique-template deletion and reference-spanning evidence."""
    del_templates = set()
    ref_templates = set()
    breakpoint_pairs = defaultdict(list)
    scan_start = max(0, core_start - scan_flank)
    scan_end = core_end + scan_flank
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, scan_start, scan_end):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_qcfail
                or read.is_duplicate
                or read.mapping_quality < 20
            ):
                continue
            template = read.query_name
            segments = sorted(_sa_segments(read, chrom, 20))
            blocks = _cigar_blocks(read)
            if continuous_core_span(blocks, core_start, core_end):
                ref_templates.add(template)
            for left, right in zip(segments, segments[1:]):
                if gap_spans_core(left[1], right[0], core_start, core_end):
                    del_templates.add(template)
                    breakpoint_pairs[template].append((left[1], right[0]))
        flank_depths = [
            mean_depth(bam, chrom, core_start - 500, core_start),
            mean_depth(bam, chrom, core_end, core_end + 500),
        ]
        flank_depth = sum(flank_depths) / len(flank_depths)
        core_depth = mean_depth(bam, chrom, core_start, core_end)
    return del_templates, ref_templates, breakpoint_pairs, flank_depth, core_depth


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chrom", default="Chr4")
    parser.add_argument("--core-start", type=int, default=CORE_START)
    parser.add_argument("--core-end", type=int, default=CORE_END)
    args = parser.parse_args()

    del_templates, ref_templates, breakpoint_pairs, flank_depth, core_depth = collect_evidence(
        args.bam, args.chrom, args.core_start, args.core_end
    )
    pair_counts = Counter(pair for pairs in breakpoint_pairs.values() for pair in pairs)
    modal = "."
    if pair_counts:
        (left, right), n = pair_counts.most_common(1)[0]
        modal = f"{left + 1}-{right}({n})"
    genotype = call_genotype(flank_depth, core_depth, len(ref_templates), len(del_templates))
    fields = [
        "sample", "chrom", "core_start", "core_end", "flank_mean_depth", "core_mean_depth",
        "core_to_flank_ratio", "ref_continuous_templates", "del_split_templates",
        "modal_del_breakpoint_1based", "GT", "evidence_class",
    ]
    row = {
        "sample": args.sample,
        "chrom": args.chrom,
        "core_start": args.core_start + 1,
        "core_end": args.core_end,
        "flank_mean_depth": f"{flank_depth:.4f}",
        "core_mean_depth": f"{core_depth:.4f}",
        "core_to_flank_ratio": f"{core_depth / flank_depth:.4f}" if flank_depth else ".",
        "ref_continuous_templates": len(ref_templates),
        "del_split_templates": len(del_templates),
        "modal_del_breakpoint_1based": modal,
        "GT": genotype,
        "evidence_class": "HIGH_CONFIDENCE" if genotype != "./." else "INSUFFICIENT_OR_CONFLICTING",
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
