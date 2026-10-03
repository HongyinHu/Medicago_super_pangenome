#!/usr/bin/env python3
import argparse
import gzip
import os
import re
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_info(info):
    d = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def svlen_ok(info, pos, min_len):
    if "SVLEN" in info:
        vals = []
        for x in str(info["SVLEN"]).split(","):
            try:
                vals.append(abs(int(float(x))))
            except ValueError:
                pass
        if vals:
            return max(vals) >= min_len
    if "END" in info:
        try:
            return abs(int(info["END"]) - int(pos)) >= min_len
        except ValueError:
            return False
    return False


def read_fai(path):
    contigs = set()
    with open(path) as fh:
        for line in fh:
            if line.strip():
                contigs.add(line.split("\t", 1)[0])
    return contigs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-vcf", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--caller", required=True)
    ap.add_argument("--ref-fai", required=True)
    ap.add_argument("--min-svlen", type=int, default=50)
    ap.add_argument("--svtypes", default="DEL,INS,DUP,INV")
    args = ap.parse_args()

    allowed = set(x.strip() for x in args.svtypes.split(",") if x.strip())
    ref_contigs = read_fai(args.ref_fai)
    header_contigs = set()
    kept = 0
    seen = 0

    os.makedirs(os.path.dirname(args.out_vcf), exist_ok=True)
    with open_text(args.in_vcf) as inp, open(args.out_vcf, "w") as out:
        for line in inp:
            if line.startswith("##contig="):
                m = re.search(r"ID=([^,>]+)", line)
                if m:
                    header_contigs.add(m.group(1))
                out.write(line)
                continue
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                bad = sorted(header_contigs - ref_contigs)
                if bad:
                    sys.stderr.write("ERROR: %s has contigs not present in %s; first mismatches: %s\n" %
                                     (args.in_vcf, args.ref_fai, ",".join(bad[:20])))
                    return 3
                out.write("##panSV_filter=PASS_only;SVLEN>=%d;SVTYPE=%s;caller=%s\n" %
                          (args.min_svlen, ",".join(sorted(allowed)), args.caller))
                cols = line.rstrip("\n").split("\t")
                if len(cols) > 9:
                    cols = cols[:9] + [args.sample]
                out.write("\t".join(cols) + "\n")
                continue
            if not line.strip():
                continue
            seen += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, _ref, _alt, _qual, flt, info_s = cols[:8]
            if chrom not in ref_contigs:
                sys.stderr.write("ERROR: record contig %s in %s is not present in %s\n" %
                                 (chrom, args.in_vcf, args.ref_fai))
                return 3
            if flt != "PASS":
                continue
            info = parse_info(info_s)
            svtype = str(info.get("SVTYPE", "")).split(",", 1)[0]
            if svtype not in allowed:
                continue
            if not svlen_ok(info, pos, args.min_svlen):
                continue
            if cols[2] == ".":
                cols[2] = "%s_%s_%d" % (args.sample, args.caller, seen)
            if len(cols) > 10:
                cols = cols[:10]
            out.write("\t".join(cols) + "\n")
            kept += 1

    sys.stderr.write("filtered\t%s\tseen=%d\tkept=%d\n" % (args.in_vcf, seen, kept))
    if kept == 0:
        sys.stderr.write("WARNING: no records kept for %s\n" % args.in_vcf)
    return 0


if __name__ == "__main__":
    sys.exit(main())
