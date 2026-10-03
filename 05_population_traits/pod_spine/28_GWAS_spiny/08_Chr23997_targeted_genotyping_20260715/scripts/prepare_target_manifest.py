#!/usr/bin/env python3
"""Build the exact 143-sample manifest in the established GWAS order."""

import argparse
import csv
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--keep", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--bam-list", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    phenotype = {}
    with open(args.phenotype, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            phenotype[row["IID"]] = row["Pod_spine_binary"]

    bams = {}
    with open(args.bam_list) as handle:
        for line in handle:
            bam = line.strip()
            if not bam:
                continue
            basename = os.path.basename(bam)
            if not basename.endswith(".dedup.bam"):
                raise ValueError(f"Unexpected BAM name: {bam}")
            sample = basename[: -len(".dedup.bam")]
            bams[sample] = bam

    ordered = []
    with open(args.keep) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2 or fields[0] != fields[1]:
                raise ValueError(f"Unexpected keep row: {line.rstrip()}")
            ordered.append(fields[1])

    missing_pheno = [sample for sample in ordered if sample not in phenotype]
    missing_bam = [sample for sample in ordered if sample not in bams]
    if missing_pheno or missing_bam:
        raise ValueError(
            f"Missing phenotype={missing_pheno}; missing BAM={missing_bam}"
        )
    if len(ordered) != 143 or len(set(ordered)) != 143:
        raise ValueError(f"Expected 143 unique samples, found {len(ordered)}")

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["task_index", "sample", "phenotype", "bam"])
        for index, sample in enumerate(ordered, 1):
            writer.writerow([index, sample, phenotype[sample], bams[sample]])


if __name__ == "__main__":
    main()
