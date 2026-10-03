#!/usr/bin/env python3
import argparse
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def read_fai(path):
    lengths = {}
    with open(path) as fh:
        for line in fh:
            if line.strip():
                parts = line.rstrip("\n").split("\t")
                lengths[parts[0]] = int(parts[1])
    return lengths


def parse_info(info):
    out = {}
    flags = []
    for item in info.split(";"):
        if not item or item == ".":
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            flags.append(item)
            out[item] = True
    return out, flags


def set_info(info, key, value):
    info[key] = value


def info_to_string(info, flags):
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


def int_or_none(value):
    try:
        return int(float(str(value).split(",", 1)[0]))
    except Exception:
        return None


def parse_breakend_alt(alt):
    m = re.search(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]", alt)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def mate_from_record(chrom, pos, alt, info):
    chr2 = info.get("CHR2")
    matepos = int_or_none(info.get("MATEPOS"))
    if chr2 is None or chr2 is True or matepos is None:
        alt_chr, alt_pos = parse_breakend_alt(alt)
        if chr2 is None or chr2 is True:
            chr2 = alt_chr
        if matepos is None:
            matepos = alt_pos
    return chr2, matepos


def record_end(pos, info, svtype):
    end = int_or_none(info.get("END"))
    if end is not None:
        return max(pos, end)
    svlen = int_or_none(info.get("SVLEN"))
    if svlen is not None and svtype not in ("INS", "TRA", "BND"):
        return max(pos, pos + abs(svlen))
    return pos


def add_window(windows, contig_lengths, var_id, chrom, center, flank, tag):
    if chrom not in contig_lengths or center is None:
        return False
    clen = contig_lengths[chrom]
    c0 = max(0, min(clen, int(center) - 1))
    if tag.endswith("left"):
        start = max(0, c0 - flank)
        end = c0
    else:
        start = min(clen, c0 + 1)
        end = min(clen, c0 + 1 + flank)
    if end <= start:
        return False
    name = "%s|%s" % (var_id, tag)
    windows.append((chrom, start, end, name))
    return True


def build_windows(record, idx, contig_lengths, flank):
    cols = record.rstrip("\n").split("\t")
    chrom = cols[0]
    pos = int(cols[1])
    rec_id = cols[2] if cols[2] != "." else "record_%d" % idx
    alt = cols[4]
    info, _flags = parse_info(cols[7])
    svtype = str(info.get("SVTYPE", ".")).split(",", 1)[0]
    end = record_end(pos, info, svtype)
    windows = []
    if svtype == "TRA" or svtype == "BND":
        ok1 = add_window(windows, contig_lengths, rec_id, chrom, pos, flank, "local_left")
        ok2 = add_window(windows, contig_lengths, rec_id, chrom, pos, flank, "local_right")
        chr2, matepos = mate_from_record(chrom, pos, alt, info)
        ok3 = add_window(windows, contig_lengths, rec_id, chr2, matepos, flank, "mate_left")
        ok4 = add_window(windows, contig_lengths, rec_id, chr2, matepos, flank, "mate_right")
        return rec_id, windows, ok1 and ok2 and ok3 and ok4
    ok1 = add_window(windows, contig_lengths, rec_id, chrom, pos, flank, "left")
    ok2 = add_window(windows, contig_lengths, rec_id, chrom, end, flank, "right")
    return rec_id, windows, ok1 and ok2


def run_bedcov(samtools, bam, bed, out_path):
    with open(out_path, "w") as out:
        subprocess.check_call([samtools, "bedcov", bed, bam], stdout=out)


def run_mosdepth(mosdepth, threads, bam, bed, prefix):
    subprocess.check_call([mosdepth, "-t", str(threads), "--by", bed, prefix, bam])
    return prefix + ".regions.bed.gz"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-vcf", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--qc-tsv", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref-fai", required=True)
    ap.add_argument("--samtools", required=True)
    ap.add_argument("--mosdepth", default="")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--tmp-prefix", required=True)
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--min-each-flank-depth", type=float, default=5.0)
    ap.add_argument("--min-mean-flank-depth", type=float, default=8.0)
    args = ap.parse_args()

    contig_lengths = read_fai(args.ref_fai)
    os.makedirs(os.path.dirname(args.out_vcf), exist_ok=True)
    os.makedirs(os.path.dirname(args.qc_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.tmp_prefix), exist_ok=True)

    header = []
    records = []
    windows = []
    rec_windows_ok = {}
    rec_order = []
    with open_text(args.in_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            if not line.strip():
                continue
            idx = len(records) + 1
            rec_id, rec_windows, ok = build_windows(line, idx, contig_lengths, args.flank)
            records.append((rec_id, line))
            rec_order.append(rec_id)
            rec_windows_ok[rec_id] = ok
            windows.extend(rec_windows)

    bed_path = args.tmp_prefix + ".flanks.bed"
    bedcov_path = args.tmp_prefix + ".flanks.bedcov.tsv"
    with open(bed_path, "w") as bed:
        for chrom, start, end, name in windows:
            bed.write("%s\t%d\t%d\t%s\n" % (chrom, start, end, name))

    depths = defaultdict(list)
    if windows:
        if args.mosdepth:
            regions_path = run_mosdepth(args.mosdepth, args.threads, args.bam, bed_path, args.tmp_prefix + ".mosdepth")
            opener = gzip.open if regions_path.endswith(".gz") else open
            with opener(regions_path, "rt", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    name = parts[3]
                    try:
                        mean_depth = float(parts[-1])
                    except ValueError:
                        continue
                    var_id = name.split("|", 1)[0]
                    depths[var_id].append(mean_depth)
        else:
            run_bedcov(args.samtools, args.bam, bed_path, bedcov_path)
            with open(bedcov_path) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    chrom, start_s, end_s, name = parts[:4]
                    try:
                        start = int(start_s)
                        end = int(end_s)
                        cov_sum = float(parts[4])
                    except ValueError:
                        continue
                    length = max(1, end - start)
                    var_id = name.split("|", 1)[0]
                    depths[var_id].append(cov_sum / length)

    pass_ids = set()
    qc = {}
    for rec_id in rec_order:
        vals = depths.get(rec_id, [])
        if not vals or not rec_windows_ok.get(rec_id, False):
            status = "NO_DEPTH_INFO"
            min_depth = 0.0
            mean_depth = 0.0
        else:
            min_depth = min(vals)
            mean_depth = sum(vals) / len(vals)
            if min_depth >= args.min_each_flank_depth and mean_depth >= args.min_mean_flank_depth:
                status = "PASS_FLANK_QC"
                pass_ids.add(rec_id)
            else:
                status = "LOW_FLANK_COVERAGE"
        qc[rec_id] = (status, min_depth, mean_depth, len(vals))

    with open(args.out_vcf, "w") as out:
        inserted = False
        for line in header:
            if line.startswith("#CHROM") and not inserted:
                out.write("##INFO=<ID=FLANK_QC,Number=1,Type=String,Description=\"Flanking read-depth QC status\">\n")
                out.write("##INFO=<ID=FLANK_MIN_DEPTH,Number=1,Type=Float,Description=\"Minimum mean depth among tested flank windows\">\n")
                out.write("##INFO=<ID=FLANK_MEAN_DEPTH,Number=1,Type=Float,Description=\"Mean depth across tested flank windows\">\n")
                out.write("##INFO=<ID=FLANK_WINDOWS,Number=1,Type=Integer,Description=\"Number of flank windows tested\">\n")
                inserted = True
            out.write(line)
        for rec_id, line in records:
            if rec_id not in pass_ids:
                continue
            cols = line.rstrip("\n").split("\t")
            info, flags = parse_info(cols[7])
            status, min_depth, mean_depth, nwin = qc[rec_id]
            set_info(info, "FLANK_QC", status)
            set_info(info, "FLANK_MIN_DEPTH", "%.3f" % min_depth)
            set_info(info, "FLANK_MEAN_DEPTH", "%.3f" % mean_depth)
            set_info(info, "FLANK_WINDOWS", str(nwin))
            cols[7] = info_to_string(info, flags)
            out.write("\t".join(cols) + "\n")

    with open(args.qc_tsv, "w") as out:
        out.write("variant_id\tstatus\tmin_flank_depth\tmean_flank_depth\tflank_windows\n")
        for rec_id in rec_order:
            status, min_depth, mean_depth, nwin = qc[rec_id]
            out.write("%s\t%s\t%.3f\t%.3f\t%d\n" % (rec_id, status, min_depth, mean_depth, nwin))

    sys.stderr.write("flank_qc\t%s\ttotal=%d\tpass=%d\tmin_each=%.3f\tmin_mean=%.3f\n" %
                     (args.in_vcf, len(records), len(pass_ids), args.min_each_flank_depth, args.min_mean_flank_depth))


if __name__ == "__main__":
    main()
