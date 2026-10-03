#!/usr/bin/env python3
"""Validate MATEID completeness and adjacency in a Smoove sites VCF."""

import argparse
import gzip


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def info_value(info, key):
    prefix = key + "="
    for item in info.split(";"):
        if item.startswith(prefix):
            return item[len(prefix) :]
    return ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("--output", required=True)
    parser.add_argument("--require-adjacent", action="store_true")
    args = parser.parse_args()

    mate_by_id = {}
    mates = []
    previous_id = ""
    adjacent_pairs = 0

    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or info_value(fields[7], "SVTYPE") != "BND":
                previous_id = fields[2] if len(fields) > 2 else ""
                continue
            record_id = fields[2]
            mate_id = info_value(fields[7], "MATEID")
            mate_by_id[record_id] = mate_id
            mates.append((record_id, mate_id))
            if mate_id == previous_id:
                adjacent_pairs += 1
            previous_id = record_id

    missing = [(record_id, mate_id) for record_id, mate_id in mates if not mate_id or mate_id not in mate_by_id]
    nonreciprocal = [
        (record_id, mate_id)
        for record_id, mate_id in mates
        if mate_id in mate_by_id and mate_by_id[mate_id] != record_id
    ]
    reciprocal = len(mates) - len(missing) - len(nonreciprocal)

    with open(args.output, "w", encoding="utf-8", newline="") as out:
        out.write("metric\tvalue\n")
        out.write(f"bnd_records\t{len(mates)}\n")
        out.write(f"unique_bnd_ids\t{len(mate_by_id)}\n")
        out.write(f"records_with_reciprocal_mate\t{reciprocal}\n")
        out.write(f"records_with_missing_mate\t{len(missing)}\n")
        out.write(f"records_with_nonreciprocal_mate\t{len(nonreciprocal)}\n")
        out.write(f"adjacent_completed_pairs\t{adjacent_pairs}\n")

    if missing:
        raise SystemExit(f"Found {len(missing)} BND records with absent MATEID partners")
    if nonreciprocal:
        raise SystemExit(f"Found {len(nonreciprocal)} BND records with non-reciprocal MATEID partners")
    if args.require_adjacent and adjacent_pairs * 2 != len(mates):
        raise SystemExit(
            f"Only {adjacent_pairs} of {len(mates) // 2} BND pairs are adjacent"
        )


if __name__ == "__main__":
    main()
