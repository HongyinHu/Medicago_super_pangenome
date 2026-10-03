#!/usr/bin/env python3
"""Genotype the Chr23997 intron-2 deletion by local two-haplotype alignment.

This script produces a conservative genotype only when reads discriminate the
intact reference haplotype from the exact deletion haplotype.  The breakpoint
coordinates follow the candidate VCF: POS=88,980,831 and END=88,981,052 on
genome_Msa Chr4, therefore the deleted sequence is 221 bp.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pysam


def build_alt_haplotype(sequence: str, deletion_start: int, deletion_end: int):
    """Remove a 0-based half-open interval and return sequence plus junction."""
    if not (0 <= deletion_start < deletion_end <= len(sequence)):
        raise ValueError("invalid deletion interval")
    return sequence[:deletion_start] + sequence[deletion_end:], deletion_start


def junction_kmers(sequence: str, junction: int, kmer_size: int):
    """Return k-mers that contain bases on both sides of a deletion junction."""
    values = []
    for start in range(max(0, junction - kmer_size + 1), min(junction, len(sequence) - kmer_size + 1)):
        end = start + kmer_size
        if start < junction < end:
            values.append(sequence[start:end].upper())
    return values


def call_genotype(flank_depth: float, ref_reads: int, alt_reads: int):
    """Return a conservative diploid genotype from allele-specific read counts."""
    if flank_depth < 8 or max(ref_reads, alt_reads) < 3:
        return "./."
    if ref_reads >= 3 and alt_reads >= 3:
        return "0/1"
    if ref_reads >= 3 and alt_reads == 0:
        return "0/0"
    if alt_reads >= 3 and ref_reads == 0:
        return "1/1"
    return "./."


def reverse_complement(sequence: str):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def parse_alignment_scores(path: Path, ref_name: str, alt_name: str, ref_inner_start: int,
                           ref_inner_end: int, alt_junction: int, min_anchor: int):
    """Summarize allele-discriminating local BWA alignments by read template."""
    per_read = defaultdict(dict)
    with pysam.AlignmentFile(path, "r") as sam:
        for read in sam:
            # Both competing haplotypes are reported as primary/secondary BWA
            # alignments.  Keep secondaries here, but retain only the best score
            # per haplotype below.
            if read.is_unmapped or read.is_supplementary:
                continue
            if read.mapping_quality < 20 or not read.has_tag("AS"):
                continue
            start = read.reference_start
            end = read.reference_end or start
            record = {
                "as": int(read.get_tag("AS")),
                "start": start,
                "end": end,
                "mapq": int(read.mapping_quality),
            }
            raw_allele = read.reference_name
            if raw_allele == ref_name or raw_allele.startswith(ref_name + "|"):
                allele = ref_name
            elif raw_allele == alt_name or raw_allele.startswith(alt_name + "|"):
                allele = alt_name
            else:
                continue
            previous = per_read[read.query_name].get(allele)
            if previous is None or record["as"] > previous["as"]:
                per_read[read.query_name][allele] = record

    ref_support = set()
    alt_support = set()
    for query_name, values in per_read.items():
        ref = values.get(ref_name)
        alt = values.get(alt_name)
        if ref is None and alt is None:
            continue
        ref_as = ref["as"] if ref else -10_000
        alt_as = alt["as"] if alt else -10_000
        if (
            alt is not None
            and alt_as - ref_as >= 10
            and alt["start"] <= alt_junction - min_anchor
            and alt["end"] >= alt_junction + min_anchor
        ):
            alt_support.add(query_name)
        if (
            ref is not None
            and ref_as - alt_as >= 10
            and ref["end"] > ref_inner_start + min_anchor
            and ref["start"] < ref_inner_end - min_anchor
        ):
            ref_support.add(query_name)
    return ref_support, alt_support


def source_flank_depth(bam_path: str, chrom: str, start0: int, end0: int, flank: int):
    values = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for left, right in ((max(0, start0 - flank), start0), (end0, end0 + flank)):
            for column in bam.pileup(chrom, left, right, truncate=True, stepper="samtools", min_mapping_quality=20):
                if left <= column.reference_pos < right:
                    values.append(column.nsegments)
    return sum(values) / len(values) if values else 0.0


def write_fastq_from_bam(bam_path: str, chrom: str, start0: int, end0: int, out_fastq: Path):
    """Extract primary high-quality reads overlapping the local haplotype window."""
    seen = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam, out_fastq.open("w") as handle:
        for read in bam.fetch(chrom, start0, end0):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_duplicate
                or read.is_qcfail
                or read.mapping_quality < 20
                or not read.query_sequence
                or read.query_name in seen
            ):
                continue
            seen.add(read.query_name)
            sequence = read.query_sequence
            qualities = read.qual or "I" * len(sequence)
            handle.write(f"@{read.query_name}\n{sequence}\n+\n{qualities}\n")
    return len(seen)


def count_exact_kmer_reads(fastq_path: Path, kmers):
    """Count templates containing an ALT-junction or REF-inner diagnostic k-mer."""
    values = {kmer.upper() for kmer in kmers}
    if not values:
        return 0
    count = 0
    with fastq_path.open() as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip().upper()
            handle.readline()
            handle.readline()
            rev = reverse_complement(sequence)
            if any(kmer in sequence or kmer in rev for kmer in values):
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--haplotypes", required=True)
    parser.add_argument("--bwa", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--tmpdir", required=True)
    parser.add_argument("--chrom", default="Chr4")
    parser.add_argument("--pos", type=int, default=88980831)
    parser.add_argument("--end", type=int, default=88981052)
    parser.add_argument("--window", type=int, default=2000)
    parser.add_argument("--kmer", type=int, default=31)
    args = parser.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(args.tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    pos0 = args.pos
    end0 = args.end
    region_start0 = max(0, pos0 - args.window)
    region_end0 = end0 + args.window

    with pysam.FastaFile(args.reference) as fasta:
        reference_haplotype = fasta.fetch(args.chrom, region_start0, region_end0).upper()
    local_del_start = pos0 - region_start0
    local_del_end = end0 - region_start0
    alternate_haplotype, alt_junction = build_alt_haplotype(
        reference_haplotype, local_del_start, local_del_end
    )

    haplotypes = Path(args.haplotypes)
    if not haplotypes.exists():
        haplotypes.parent.mkdir(parents=True, exist_ok=True)
        with haplotypes.open("w") as handle:
            handle.write(f">Chr23997_REF\n{reference_haplotype}\n")
            handle.write(f">Chr23997_DEL221\n{alternate_haplotype}\n")
        subprocess.run([args.bwa, "index", str(haplotypes)], check=True)

    fastq = tmpdir / f"{args.sample}.local.fq"
    sam = tmpdir / f"{args.sample}.haplotypes.sam"
    n_reads = write_fastq_from_bam(args.bam, args.chrom, region_start0, region_end0, fastq)
    with sam.open("w") as handle:
        subprocess.run(
            [args.bwa, "mem", "-a", "-T", "0", "-k", "19", str(haplotypes), str(fastq)],
            check=True,
            stdout=handle,
            stderr=subprocess.DEVNULL,
        )

    ref_support, alt_support = parse_alignment_scores(
        sam,
        "Chr23997_REF",
        "Chr23997_DEL221",
        local_del_start,
        local_del_end,
        alt_junction,
        20,
    )
    alt_kmers = junction_kmers(alternate_haplotype, alt_junction, args.kmer)
    inner = reference_haplotype[local_del_start:local_del_end]
    ref_kmers = [inner[start:start + args.kmer] for start in range(0, max(0, len(inner) - args.kmer + 1), 10)]
    alt_kmer_reads = count_exact_kmer_reads(fastq, alt_kmers)
    ref_kmer_reads = count_exact_kmer_reads(fastq, ref_kmers)
    flank_depth = source_flank_depth(args.bam, args.chrom, pos0, end0, 500)
    genotype = call_genotype(flank_depth, len(ref_support), len(alt_support))

    fields = [
        "sample", "chrom", "pos", "end", "svlen", "flank_mean_depth", "local_templates",
        "ref_allele_reads", "alt_allele_reads", "ref_inner_kmer_templates", "alt_junction_kmer_templates",
        "GT", "evidence_class",
    ]
    row = {
        "sample": args.sample,
        "chrom": args.chrom,
        "pos": args.pos,
        "end": args.end,
        "svlen": -(args.end - args.pos),
        "flank_mean_depth": f"{flank_depth:.4f}",
        "local_templates": n_reads,
        "ref_allele_reads": len(ref_support),
        "alt_allele_reads": len(alt_support),
        "ref_inner_kmer_templates": ref_kmer_reads,
        "alt_junction_kmer_templates": alt_kmer_reads,
        "GT": genotype,
        "evidence_class": "HIGH_CONFIDENCE" if genotype != "./." else "INSUFFICIENT_OR_CONFLICTING",
    }
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
