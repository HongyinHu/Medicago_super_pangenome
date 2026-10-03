#!/usr/bin/env python3
"""Validate the genome_Msa paper-style unified panSV outputs."""

import csv
from collections import Counter
from pathlib import Path


RUN_DIR = Path("path/to/project/N_3.call_SV/06.integrated_panSV_Msa_paper_style_unified")
PREFIX = "paper_style_unified_panSV.Msa_ref"
FINAL_SAMPLES = [
    "genome_395",
    "genome_410",
    "genome_436",
    "genome_454",
    "genome_457",
    "genome_461",
    "genome_468",
    "genome_472",
    "genome_474",
    "genome_482",
    "genome_M22",
    "genome_M46",
    "genome_Mar",
    "genome_Mpo",
    "genome_Mru",
    "genome_Msa",
    "genome_R108",
    "genome_ZM4",
]


def count_lines(path):
    n = 0
    with path.open("rb") as handle:
        for _ in handle:
            n += 1
    return n


def main():
    results = RUN_DIR / "results"
    summary = RUN_DIR / "summary"
    events = results / (PREFIX + ".events.tsv")
    pav = results / (PREFIX + ".PAV.matrix.tsv")
    stats = results / (PREFIX + ".stats.tsv")
    evidence = results / (PREFIX + ".svgap_readbased.overlap_evidence.tsv")
    counts_file = summary / (PREFIX + ".counts.tsv")
    dropped = summary / (PREFIX + ".dropped_excluded_samples.tsv")
    done = summary / (PREFIX + ".done")

    paths = [events, pav, stats, evidence, counts_file, dropped, done]
    print("FILES")
    for path in paths:
        size_mb = path.stat().st_size / 1024 / 1024
        print("%s\t%.2f MB\t%s lines" % (path, size_mb, count_lines(path) if path.suffix == ".tsv" else "NA"))

    source_counts = Counter()
    confidence_counts = Counter()
    svtype_counts = Counter()
    bad_labels = Counter()
    msa_prefix = 0
    svgap_genome474_present = 0

    with events.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            source_counts[row["source_methods"]] += 1
            confidence_counts[row["final_confidence"]] += 1
            svtype_counts[row["SVTYPE"]] += 1
            if row["CHROM"].startswith("Msa."):
                msa_prefix += 1
            blob = "\t".join(row.get(k, "") for k in row.keys())
            for label in ("genome_A17", "A17", "genome_G474b", "G474b", "genome_G474a", "G474a"):
                if label in blob:
                    bad_labels[label] += 1
            if row["source_methods"] == "SVGAP" and "genome_474" in row.get("svgap_support_samples", "").split(","):
                svgap_genome474_present += 1

    with pav.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        missing_samples = [sample for sample in FINAL_SAMPLES if sample not in header]
        unexpected_genome_columns = [
            col for col in header
            if col.startswith("genome_") and col not in FINAL_SAMPLES
        ]
        pav_rows = 0
        pav_genome474_ones = 0
        pav_svgap_genome474_ones = 0
        pav_bad_chrom = 0
        for row in reader:
            pav_rows += 1
            if row["CHROM"].startswith("Msa."):
                pav_bad_chrom += 1
            if row["genome_474"] == "1":
                pav_genome474_ones += 1
                if row["source_methods"] == "SVGAP":
                    pav_svgap_genome474_ones += 1

    print("SUMMARY")
    print("event_rows\t%s" % sum(source_counts.values()))
    print("pav_rows\t%s" % pav_rows)
    print("source_methods\t%s" % dict(source_counts))
    print("final_confidence\t%s" % dict(confidence_counts))
    print("SVTYPE\t%s" % dict(svtype_counts))
    print("missing_sample_columns\t%s" % ",".join(missing_samples))
    print("unexpected_genome_columns\t%s" % ",".join(unexpected_genome_columns))
    print("bad_label_hits\t%s" % dict(bad_labels))
    print("events_with_Msa_prefix_chrom\t%s" % msa_prefix)
    print("pav_rows_with_Msa_prefix_chrom\t%s" % pav_bad_chrom)
    print("SVGAP_only_events_with_genome_474_support\t%s" % svgap_genome474_present)
    print("PAV_genome_474_presence_rows\t%s" % pav_genome474_ones)
    print("PAV_SVGAP_only_genome_474_presence_rows\t%s" % pav_svgap_genome474_ones)


if __name__ == "__main__":
    main()
