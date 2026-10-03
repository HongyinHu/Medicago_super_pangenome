#!/usr/bin/env python3
"""Filter Smoove/SVtools cohort sites without performing any SV clustering."""

import argparse
import gzip
import re


PRIMARY_CONTIGS = {f"Chr{i}" for i in range(1, 9)}
ALLOWED_TYPES = {"DEL", "DUP", "INV", "BND"}


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")


def parse_info(text):
    parsed = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = None
    return parsed


def discovery_samples(sname):
    samples = set()
    if not sname:
        return samples
    for token in re.split(r"[,|]", sname):
        token = token.strip()
        if token:
            samples.add(token.split(":", 1)[0])
    return samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-discovery-samples", type=int, default=2)
    parser.add_argument("--max-svlen", type=int, default=1_000_000)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    counts = {"input": 0, "kept": 0, "non_primary_contig": 0, "svtype": 0, "svlen": 0, "discovery_samples": 0}
    type_counts = {}

    with open_text(args.input) as source, open(args.output, "w", encoding="utf-8", newline="") as out:
        added_header = False
        for line in source:
            if line.startswith("##INFO=<ID=PRPOS") or line.startswith("##INFO=<ID=PREND"):
                continue
            if line.startswith("#CHROM") and not added_header:
                out.write('##INFO=<ID=DISCOVERY_NSAMP,Number=1,Type=Integer,Description="Number of unique samples contributing discovery evidence before cohort genotyping">\n')
                added_header = True
            if line.startswith("#"):
                out.write(line)
                continue

            counts["input"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            if fields[0] not in PRIMARY_CONTIGS:
                counts["non_primary_contig"] += 1
                continue

            info = parse_info(fields[7])
            svtype = info.get("SVTYPE", "")
            if svtype not in ALLOWED_TYPES:
                counts["svtype"] += 1
                continue

            if svtype != "BND":
                try:
                    svlen = abs(int(float(info.get("SVLEN", "0").split(",")[0])))
                except ValueError:
                    svlen = 0
                if svlen < 50 or svlen > args.max_svlen:
                    counts["svlen"] += 1
                    continue

            nsamp = len(discovery_samples(info.get("SNAME", "")))
            if nsamp < args.min_discovery_samples:
                counts["discovery_samples"] += 1
                continue

            cleaned = []
            for item in fields[7].split(";"):
                if item.startswith("PRPOS=") or item.startswith("PREND=") or item.startswith("DISCOVERY_NSAMP="):
                    continue
                cleaned.append(item)
            cleaned.append(f"DISCOVERY_NSAMP={nsamp}")
            fields[7] = ";".join(cleaned)
            out.write("\t".join(fields) + "\n")
            counts["kept"] += 1
            type_counts[svtype] = type_counts.get(svtype, 0) + 1

    with open(args.summary, "w", encoding="utf-8", newline="") as out:
        out.write("metric\tvalue\n")
        for key, value in counts.items():
            out.write(f"{key}\t{value}\n")
        for svtype in sorted(type_counts):
            out.write(f"kept_{svtype}\t{type_counts[svtype]}\n")


if __name__ == "__main__":
    main()
