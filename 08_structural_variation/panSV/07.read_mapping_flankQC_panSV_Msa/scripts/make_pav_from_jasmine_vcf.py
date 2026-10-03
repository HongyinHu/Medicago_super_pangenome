#!/usr/bin/env python3
import argparse
import csv
import gzip
import re


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_info(info):
    out = {}
    for item in info.split(";"):
        if not item or item == ".":
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = True
    return out


def int_or_dot(value):
    try:
        return str(int(float(str(value).split(",", 1)[0])))
    except Exception:
        return "."


def parse_breakend_alt(alt):
    m = re.search(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]", alt)
    if not m:
        return ".", "."
    return m.group(1), m.group(2)


def read_samples(path):
    samples = []
    with open(path) as fh:
        for line in fh:
            value = line.strip()
            if value:
                samples.append(value)
    return samples


def supp_vec(info, samples):
    n = len(samples)
    # For species-level PAV, use only a vector whose length exactly matches the
    # species list. SUPP_VEC_EXT may be longer after nested/caller-level merges
    # and must not be truncated to species columns.
    for key in ("SUPP_VEC", "SUPP_VEC_EXT"):
        value = info.get(key)
        if value and value is not True:
            bits = str(value).strip()
            if len(bits) == n:
                return bits
    return None


def gt_presence(format_s, sample_s):
    if not format_s or not sample_s:
        return "0"
    fmt = format_s.split(":")
    vals = sample_s.split(":")
    if "GT" not in fmt:
        return "0"
    gt = vals[fmt.index("GT")]
    if gt in ("0/1", "1/0", "1/1", "0|1", "1|0", "1|1"):
        return "1"
    if gt in ("./.", ".|."):
        return "NA"
    return "0"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--species-list", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--stats", required=True)
    args = ap.parse_args()
    samples = read_samples(args.species_list)
    counts = {}

    with open_text(args.vcf) as inp, open(args.matrix, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["SV_ID", "CHROM", "POS", "END", "SVTYPE", "SVLEN", "CHR2", "MATEPOS"] + samples)
        header_samples = []
        for line in inp:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                header_samples = parts[9:] if len(parts) > 9 else []
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, rec_id, _ref, alt, _qual, _flt, info_s = parts[:8]
            info = parse_info(info_s)
            svtype = str(info.get("SVTYPE", ".")).split(",", 1)[0]
            if svtype == "BND":
                svtype = "TRA"
            end = int_or_dot(info.get("END", pos))
            svlen = int_or_dot(info.get("SVLEN", "."))
            chr2 = info.get("CHR2", ".")
            matepos = int_or_dot(info.get("MATEPOS", "."))
            if svtype == "TRA" and (chr2 == "." or matepos == "."):
                alt_chr, alt_pos = parse_breakend_alt(alt)
                chr2 = alt_chr
                matepos = alt_pos
            bits = supp_vec(info, samples)
            calls = []
            if bits is not None:
                calls = ["1" if bit == "1" else "0" for bit in bits]
            elif header_samples and len(parts) >= 10:
                sample_to_call = {}
                fmt = parts[8]
                for sample, value in zip(header_samples, parts[9:]):
                    sample_to_call[sample] = gt_presence(fmt, value)
                calls = [sample_to_call.get(sample, "0") for sample in samples]
            else:
                calls = ["0"] * len(samples)
            counts[svtype] = counts.get(svtype, 0) + 1
            writer.writerow([rec_id, chrom, pos, end, svtype, svlen, chr2, matepos] + calls)

    with open(args.stats, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["section", "key", "count"])
        writer.writerow(["input", "samples", len(samples)])
        for key in sorted(counts):
            writer.writerow(["SVTYPE", key, counts[key]])
        writer.writerow(["rows", "events", sum(counts.values())])


if __name__ == "__main__":
    main()
