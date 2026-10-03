#!/usr/bin/env python3
import argparse
import gzip
import os
from collections import Counter, defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_info(info):
    d = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def gt_state(format_s, sample_s):
    keys = format_s.split(":")
    vals = sample_s.split(":")
    d = {k: vals[i] if i < len(vals) else "." for i, k in enumerate(keys)}
    gt = d.get("GT", "./.")
    if gt in ("./.", ".|.", "."):
        return "NA"
    alleles = gt.replace("|", "/").split("/")
    if all(a == "0" for a in alleles if a != "."):
        return "0"
    if any(a not in ("0", ".") for a in alleles):
        return "1"
    return "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--stats", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.matrix), exist_ok=True)
    samples = []
    rows = []
    svtype_by_id = {}
    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            info = parse_info(c[7])
            sv_id = c[2]
            svtype = info.get("SVTYPE", "NA").split(",", 1)[0]
            end = info.get("END", ".")
            svlen = info.get("SVLEN", ".")
            states = [gt_state(c[8], x) for x in c[9:]]
            rows.append([sv_id, c[0], c[1], end, svtype, svlen] + states)
            svtype_by_id[sv_id] = svtype

    with open(args.matrix, "w") as out:
        out.write("SV_ID\tCHROM\tPOS\tEND\tSVTYPE\tSVLEN\t%s\n" % "\t".join(samples))
        for row in rows:
            out.write("\t".join(row) + "\n")

    stats = {s: Counter() for s in samples}
    type_counts = {s: Counter() for s in samples}
    for row in rows:
        sv_id = row[0]
        svtype = svtype_by_id.get(sv_id, "NA")
        vals = row[6:]
        present = [samples[i] for i, v in enumerate(vals) if v == "1"]
        for i, s in enumerate(samples):
            v = vals[i]
            if v == "1":
                stats[s]["present"] += 1
                type_counts[s][svtype] += 1
            elif v == "0":
                stats[s]["absent"] += 1
            else:
                stats[s]["missing"] += 1
        if len(present) == 1:
            stats[present[0]]["singleton_present"] += 1
    all_types = sorted(set(t for c in type_counts.values() for t in c))
    with open(args.stats, "w") as out:
        out.write("sample\tpresent\tabsent\tmissing\tmissing_rate\tsingleton_present\t%s\n" %
                  "\t".join("SVTYPE_" + t for t in all_types))
        total = float(len(rows)) if rows else 1.0
        for s in samples:
            miss = stats[s]["missing"]
            out.write("%s\t%d\t%d\t%d\t%.6f\t%d\t%s\n" % (
                s, stats[s]["present"], stats[s]["absent"], miss, miss / total,
                stats[s]["singleton_present"],
                "\t".join(str(type_counts[s][t]) for t in all_types)))


if __name__ == "__main__":
    main()
