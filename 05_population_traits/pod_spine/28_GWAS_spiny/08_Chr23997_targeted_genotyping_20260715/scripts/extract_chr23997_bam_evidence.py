#!/usr/bin/env python3
"""Extract independent short-read evidence for the validated Chr23997 deletion."""

import argparse
import csv
import math
import re
from collections import defaultdict

import pysam


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def depth_vector(bam, chrom, start, end, min_mapq):
    depths = [0] * max(0, end - start)
    for col in bam.pileup(
        chrom,
        start,
        end,
        truncate=True,
        stepper="samtools",
        min_base_quality=0,
        min_mapping_quality=min_mapq,
        max_depth=1000000,
    ):
        if start <= col.reference_pos < end:
            depths[col.reference_pos - start] = col.nsegments
    return depths


def alignment_blocks_and_deletions(read):
    ref = read.reference_start
    deletions = []
    if not read.cigartuples:
        return deletions
    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X
            ref += length
        elif op == 2:  # D
            deletions.append((ref, ref + length, length))
            ref += length
        elif op == 3:  # N
            ref += length
        elif op in (1, 4, 5, 6):
            continue
    return deletions


def median(values):
    if not values:
        return float("nan")
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[middle])
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def cigar_reference_length(cigar):
    return sum(
        int(length)
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if op in "MDN=X"
    )


def mate_key(read):
    if read.is_read1:
        return (read.query_name, 1)
    if read.is_read2:
        return (read.query_name, 2)
    return (read.query_name, 0)


def spans_breakpoint(read, breakpoint0, min_anchor=20):
    if read.reference_start > breakpoint0 - min_anchor:
        return False
    if read.reference_end is None or read.reference_end < breakpoint0 + min_anchor:
        return False
    for d_start, d_end, d_len in alignment_blocks_and_deletions(read):
        if d_len >= 20 and d_start - 10 <= breakpoint0 <= d_end + 10:
            return False
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chrom", default="Chr4")
    parser.add_argument("--pos", type=int, default=88980831, help="1-based VCF POS anchor")
    parser.add_argument("--end", type=int, default=88981052, help="1-based inclusive VCF END")
    parser.add_argument("--flank", type=int, default=500)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--breakpoint-tolerance", type=int, default=100)
    args = parser.parse_args()

    # The retained anchor is at POS. Deleted sequence is POS+1..END (1-based),
    # represented by [POS, END) in pysam's 0-based half-open coordinates.
    del_start0 = args.pos
    del_end0 = args.end
    left_start = max(0, del_start0 - args.flank)
    right_end = del_end0 + args.flank

    bam = pysam.AlignmentFile(args.bam, "rb")
    left_depths = depth_vector(bam, args.chrom, left_start, del_start0, args.min_mapq)
    target_depths = depth_vector(bam, args.chrom, del_start0, del_end0, args.min_mapq)
    right_depths = depth_vector(bam, args.chrom, del_end0, right_end, args.min_mapq)

    read_sets = defaultdict(set)
    segments_by_mate = defaultdict(set)
    pair_spans = {}
    tlen_values = []
    fetch_start = max(0, left_start - 1000)
    fetch_end = right_end + 1000
    for read in bam.fetch(args.chrom, fetch_start, fetch_end):
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_qcfail
            or read.is_duplicate
            or read.mapping_quality < args.min_mapq
        ):
            continue
        qname = read.query_name
        if (
            not read.is_supplementary
            and read.is_read1
            and read.is_proper_pair
            and read.template_length
        ):
            tlen_values.append(abs(read.template_length))

        key = mate_key(read)
        if read.reference_end is not None:
            segments_by_mate[key].add(
                (read.reference_start, read.reference_end, int(read.is_supplementary))
            )
        if read.has_tag("SA"):
            for item in read.get_tag("SA").rstrip(";").split(";"):
                fields = item.split(",")
                if len(fields) < 6 or fields[0] != args.chrom:
                    continue
                sa_start = int(fields[1]) - 1
                sa_end = sa_start + cigar_reference_length(fields[3])
                segments_by_mate[key].add((sa_start, sa_end, 1))

        if (
            not read.is_supplementary
            and read.is_read1
            and read.is_paired
            and not read.mate_is_unmapped
            and read.next_reference_id == read.reference_id
            and read.template_length
        ):
            template_start = min(read.reference_start, read.next_reference_start)
            pair_spans[qname] = (
                abs(read.template_length),
                template_start,
                template_start + abs(read.template_length),
            )

        for d_start, d_end, d_len in alignment_blocks_and_deletions(read):
            if 150 <= d_len <= 300:
                read_sets["large_del_any"].add(qname)
                if (
                    abs(d_start - del_start0) <= args.breakpoint_tolerance
                    and abs(d_end - del_end0) <= args.breakpoint_tolerance
                ):
                    read_sets["target_del"].add(qname)

        if spans_breakpoint(read, del_start0):
            read_sets["left_ref"].add(qname)
        if spans_breakpoint(read, del_end0):
            read_sets["right_ref"].add(qname)

        if read.cigartuples:
            if read.cigartuples[0][0] in (4, 5) and read.cigartuples[0][1] >= 15:
                if abs(read.reference_start - del_end0) <= args.breakpoint_tolerance:
                    read_sets["right_clip"].add(qname)
            if read.cigartuples[-1][0] in (4, 5) and read.cigartuples[-1][1] >= 15:
                if read.reference_end is not None and abs(read.reference_end - del_start0) <= args.breakpoint_tolerance:
                    read_sets["left_clip"].add(qname)

    for key, segments in segments_by_mate.items():
        ordered = list(segments)
        for left_segment in ordered:
            if abs(left_segment[1] - del_start0) > args.breakpoint_tolerance:
                continue
            for right_segment in ordered:
                if left_segment == right_segment:
                    continue
                if abs(right_segment[0] - del_end0) > args.breakpoint_tolerance:
                    continue
                if left_segment[2] or right_segment[2]:
                    read_sets["target_split"].add(key[0])

    median_tlen = median(tlen_values)
    mad_tlen = median([abs(value - median_tlen) for value in tlen_values]) if tlen_values else float("nan")
    if tlen_values:
        robust_sigma = 1.4826 * mad_tlen
        discordant_threshold = median_tlen + max(100.0, 4.0 * robust_sigma)
        for qname, (abs_tlen, template_start, template_end) in pair_spans.items():
            if (
                abs_tlen >= discordant_threshold
                and template_start <= del_start0 - 20
                and template_end >= del_end0 + 20
            ):
                read_sets["target_discordant_pair"].add(qname)
    else:
        discordant_threshold = float("nan")

    bam.close()

    left_mean = mean(left_depths)
    target_mean = mean(target_depths)
    right_mean = mean(right_depths)
    flank_mean = mean([left_mean, right_mean])
    depth_ratio = target_mean / flank_mean if flank_mean and not math.isnan(flank_mean) else float("nan")
    covered_fraction = (
        sum(1 for value in target_depths if value >= 5) / len(target_depths)
        if target_depths
        else float("nan")
    )
    direct_alt_reads = read_sets["target_del"] | read_sets["target_split"]
    clip_reads = read_sets["left_clip"] | read_sets["right_clip"]

    if flank_mean < 5:
        evidence_call = "LOW_COVERAGE"
    elif len(direct_alt_reads) >= 2:
        evidence_call = "DEL_DIRECT_SUPPORT"
    elif len(direct_alt_reads) >= 1 and (
        len(read_sets["target_discordant_pair"]) >= 2 or len(clip_reads) >= 2
    ):
        evidence_call = "DEL_MULTIMODAL_SUPPORT"
    elif (
        len(read_sets["target_discordant_pair"]) >= 3
        and len(clip_reads) >= 2
        and depth_ratio <= 0.75
    ):
        evidence_call = "DEL_PAIR_DEPTH_CLIP_SUPPORT"
    elif (
        depth_ratio >= 0.80
        and len(read_sets["target_del"]) == 0
        and len(read_sets["left_ref"]) >= 3
        and len(read_sets["right_ref"]) >= 3
    ):
        evidence_call = "INTACT_SUPPORT"
    else:
        evidence_call = "AMBIGUOUS"

    row = {
        "sample": args.sample,
        "bam": args.bam,
        "chrom": args.chrom,
        "vcf_pos": args.pos,
        "vcf_end": args.end,
        "svlen": -(args.end - args.pos),
        "left_flank_mean_depth": f"{left_mean:.6f}",
        "target_mean_depth": f"{target_mean:.6f}",
        "right_flank_mean_depth": f"{right_mean:.6f}",
        "flank_mean_depth": f"{flank_mean:.6f}",
        "target_to_flank_depth_ratio": f"{depth_ratio:.6f}",
        "target_fraction_depth_ge5": f"{covered_fraction:.6f}",
        "target_del_cigar_reads": len(read_sets["target_del"]),
        "target_split_alignment_reads": len(read_sets["target_split"]),
        "target_direct_alt_reads": len(direct_alt_reads),
        "target_discordant_pair_reads": len(read_sets["target_discordant_pair"]),
        "large_del_cigar_reads_any": len(read_sets["large_del_any"]),
        "left_breakpoint_ref_reads": len(read_sets["left_ref"]),
        "right_breakpoint_ref_reads": len(read_sets["right_ref"]),
        "left_breakpoint_clip_reads": len(read_sets["left_clip"]),
        "right_breakpoint_clip_reads": len(read_sets["right_clip"]),
        "target_breakpoint_clip_reads": len(clip_reads),
        "local_proper_pair_count": len(tlen_values),
        "local_proper_pair_median_abs_tlen": f"{median_tlen:.3f}",
        "local_proper_pair_mad_abs_tlen": f"{mad_tlen:.3f}",
        "discordant_pair_abs_tlen_threshold": f"{discordant_threshold:.3f}",
        "evidence_call": evidence_call,
    }
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
