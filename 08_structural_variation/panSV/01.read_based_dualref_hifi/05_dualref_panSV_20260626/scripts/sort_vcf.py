#!/usr/bin/env python3
import argparse
import gzip
import os


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def read_order(fai):
    order = {}
    with open(fai) as fh:
        for i, line in enumerate(fh):
            if line.strip():
                order[line.split("\t", 1)[0]] = i
    return order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-vcf", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--fai", required=True)
    args = ap.parse_args()
    order = read_order(args.fai)
    header = []
    records = []
    with open_text(args.in_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
            elif line.strip():
                parts = line.split("\t", 2)
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                except ValueError:
                    pos = 0
                records.append((order.get(chrom, 10**12), chrom, pos, line))
    records.sort(key=lambda x: (x[0], x[2], x[1], x[3]))
    os.makedirs(os.path.dirname(args.out_vcf), exist_ok=True)
    with open(args.out_vcf, "w") as out:
        for line in header:
            out.write(line)
        for rec in records:
            out.write(rec[3])


if __name__ == "__main__":
    main()
