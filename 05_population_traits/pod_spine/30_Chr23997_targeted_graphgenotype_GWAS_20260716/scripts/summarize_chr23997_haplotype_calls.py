#!/usr/bin/env python3
"""Merge per-sample Chr23997 local haplotype calls with phenotype metadata."""

import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--per-sample-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    outputs = []
    for record in read_tsv(args.manifest):
        path = Path(args.per_sample_dir) / f"{record['sample']}.tsv"
        if not path.exists():
            continue
        rows = read_tsv(path)
        if len(rows) != 1:
            raise ValueError(f"expected one row: {path}")
        row = dict(record)
        row.update(rows[0])
        outputs.append(row)
    fields = [
        "task_index", "sample", "phenotype", "bam", "chrom", "pos", "end", "svlen",
        "flank_mean_depth", "local_templates", "ref_allele_reads", "alt_allele_reads",
        "ref_inner_kmer_templates", "alt_junction_kmer_templates", "GT", "evidence_class",
    ]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(outputs)


if __name__ == "__main__":
    main()
