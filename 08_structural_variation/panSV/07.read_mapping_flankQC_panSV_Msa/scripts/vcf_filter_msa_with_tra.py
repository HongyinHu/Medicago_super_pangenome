#!/usr/bin/env python3
import argparse
import gzip
import os
import re
import sys


TRA_TYPES = set(["BND", "TRA", "TRL", "TRANSLOCATION", "CTX"])


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_info(info):
    out = {}
    flags = []
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            flags.append(item)
            out[item] = True
    return out, flags


def info_string(info, flags):
    keys = []
    seen = set()
    for key in flags:
        if key in info and info[key] is True and key not in seen:
            keys.append(key)
            seen.add(key)
    for key in sorted(info):
        if info[key] is True:
            if key not in seen:
                keys.append(key)
                seen.add(key)
        else:
            keys.append("%s=%s" % (key, info[key]))
    return ";".join(keys) if keys else "."


def read_fai(path):
    lengths = {}
    with open(path) as fh:
        for line in fh:
            if line.strip():
                parts = line.rstrip("\n").split("\t")
                lengths[parts[0]] = int(parts[1])
    return lengths


def parse_breakend_alt(alt):
    m = re.search(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]", alt)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def mate_from_record(chrom, pos, alt, info):
    chr2 = None
    matepos = None
    for key in ("CHR2", "CHR", "MATECHROM", "MATE_CHR", "CT"):
        value = info.get(key)
        if value and value is not True:
            if ":" in str(value):
                maybe_chr, maybe_pos = str(value).split(":", 1)
                chr2 = maybe_chr
                try:
                    matepos = int(float(maybe_pos))
                except ValueError:
                    pass
            elif str(value) not in ("3to5", "5to3", "3to3", "5to5"):
                chr2 = str(value)
            break
    for key in ("MATEPOS", "POS2", "END2"):
        value = info.get(key)
        if value and value is not True:
            try:
                matepos = int(float(str(value).split(",", 1)[0]))
            except ValueError:
                pass
            break
    if chr2 is None or matepos is None:
        a_chr, a_pos = parse_breakend_alt(alt)
        if chr2 is None:
            chr2 = a_chr
        if matepos is None:
            matepos = a_pos
    if chr2 is None:
        chr2 = chrom
    if matepos is None:
        matepos = pos
    return chr2, matepos


def svlen_ok(info, pos, min_len):
    if "SVLEN" in info:
        vals = []
        for value in str(info["SVLEN"]).split(","):
            try:
                vals.append(abs(int(float(value))))
            except ValueError:
                pass
        if vals:
            return max(vals) >= min_len
    if "END" in info:
        try:
            return abs(int(float(str(info["END"]).split(",", 1)[0])) - int(pos)) >= min_len
        except ValueError:
            return False
    return False


def is_pass_filter(value):
    return value in ("PASS", ".")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-vcf", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--caller", required=True)
    ap.add_argument("--ref-fai", required=True)
    ap.add_argument("--min-svlen", type=int, default=50)
    ap.add_argument("--svtypes", default="DEL,INS,DUP,INV,TRA")
    ap.add_argument("--interchrom-tra-only", action="store_true")
    args = ap.parse_args()

    allowed = set(x.strip() for x in args.svtypes.split(",") if x.strip())
    ref_lengths = read_fai(args.ref_fai)
    ref_contigs = set(ref_lengths)
    header_contigs = set()
    kept = 0
    seen = 0
    dropped = {}

    os.makedirs(os.path.dirname(args.out_vcf), exist_ok=True)
    with open_text(args.in_vcf) as inp, open(args.out_vcf, "w") as out:
        wrote_info = False
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
                if not wrote_info:
                    out.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Mate chromosome for TRA/BND retained by panSV filter\">\n")
                    out.write("##INFO=<ID=MATEPOS,Number=1,Type=Integer,Description=\"Mate position for TRA/BND retained by panSV filter\">\n")
                    out.write("##panSV_filter=PASS_or_dot;SVLEN>=%d_for_interval_SV;interchrom_TRA=%s;SVTYPE=%s;caller=%s\n" %
                              (args.min_svlen, str(bool(args.interchrom_tra_only)), ",".join(sorted(allowed)), args.caller))
                    wrote_info = True
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
                dropped["malformed"] = dropped.get("malformed", 0) + 1
                continue
            chrom, pos_s, rec_id, ref, alt, qual, flt, info_s = cols[:8]
            if chrom not in ref_contigs:
                sys.stderr.write("ERROR: record contig %s in %s is not present in %s\n" %
                                 (chrom, args.in_vcf, args.ref_fai))
                return 3
            if not is_pass_filter(flt):
                dropped["nonPASS"] = dropped.get("nonPASS", 0) + 1
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                dropped["bad_pos"] = dropped.get("bad_pos", 0) + 1
                continue
            info, flags = parse_info(info_s)
            svtype = str(info.get("SVTYPE", "")).split(",", 1)[0]
            is_tra = svtype in TRA_TYPES
            if is_tra:
                if "TRA" not in allowed:
                    dropped["type"] = dropped.get("type", 0) + 1
                    continue
                chr2, matepos = mate_from_record(chrom, pos, alt, info)
                if chr2 not in ref_contigs:
                    dropped["bad_mate_contig"] = dropped.get("bad_mate_contig", 0) + 1
                    continue
                if args.interchrom_tra_only and chr2 == chrom:
                    dropped["intra_bnd"] = dropped.get("intra_bnd", 0) + 1
                    continue
                info["SVTYPE"] = "TRA"
                info["CHR2"] = chr2
                info["MATEPOS"] = str(matepos)
                if "END" not in info or info["END"] is True:
                    info["END"] = str(pos)
            else:
                if svtype not in allowed:
                    dropped["type"] = dropped.get("type", 0) + 1
                    continue
                if not svlen_ok(info, pos, args.min_svlen):
                    dropped["small"] = dropped.get("small", 0) + 1
                    continue
            if cols[2] == ".":
                cols[2] = "%s_%s_%d" % (args.sample, args.caller, seen)
            if len(cols) > 10:
                cols = cols[:10]
            cols[7] = info_string(info, flags)
            out.write("\t".join(cols) + "\n")
            kept += 1

    sys.stderr.write("filtered\t%s\tseen=%d\tkept=%d\tdropped=%s\n" %
                     (args.in_vcf, seen, kept, ",".join("%s:%s" % (k, dropped[k]) for k in sorted(dropped))))
    if kept == 0:
        sys.stderr.write("WARNING: no records kept for %s\n" % args.in_vcf)
    return 0


if __name__ == "__main__":
    sys.exit(main())
