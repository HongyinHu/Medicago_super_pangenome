#!/usr/bin/env python3
"""Remove BND records whose MATEID partner was filtered out, preserving order."""

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
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    bnd_ids = set()
    with open_text(args.input) as source:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8 and info_value(fields[7], "SVTYPE") == "BND":
                bnd_ids.add(fields[2])

    input_records = 0
    kept_records = 0
    removed_orphan_bnd = 0
    with open_text(args.input) as source, open(args.output, "w", encoding="utf-8", newline="") as out:
        for line in source:
            if line.startswith("#"):
                out.write(line)
                continue
            input_records += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 8 and info_value(fields[7], "SVTYPE") == "BND":
                mate_id = info_value(fields[7], "MATEID")
                if not mate_id or mate_id not in bnd_ids:
                    removed_orphan_bnd += 1
                    continue
            out.write(line)
            kept_records += 1

    with open(args.summary, "w", encoding="utf-8", newline="") as out:
        out.write("metric\tvalue\n")
        out.write(f"input_records\t{input_records}\n")
        out.write(f"kept_records\t{kept_records}\n")
        out.write(f"removed_orphan_bnd\t{removed_orphan_bnd}\n")


if __name__ == "__main__":
    main()
