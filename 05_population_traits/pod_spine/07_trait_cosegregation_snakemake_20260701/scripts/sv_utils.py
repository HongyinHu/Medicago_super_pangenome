#!/usr/bin/env python3
import gzip
import re


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def parse_info(info):
    result = {}
    if info in ("", "."):
        return result
    for field in info.split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
        else:
            result[field] = True
    return result


def as_int(value, default=None):
    if value is None:
        return default
    try:
        return int(float(str(value).split(",")[0]))
    except ValueError:
        return default


def infer_svtype(ref, alt, info):
    svtype = info.get("SVTYPE")
    if svtype:
        return str(svtype).split(",")[0]
    m = re.search(r"<([^>]+)>", alt)
    if m:
        return m.group(1).split(":")[0]
    if len(ref) > len(alt):
        return "DEL"
    if len(alt) > len(ref):
        return "INS"
    return "UNK"


def infer_end(pos, ref, alt, info):
    end = as_int(info.get("END"))
    if end:
        return end
    svlen = abs(infer_svlen(ref, alt, info))
    if infer_svtype(ref, alt, info) == "INS":
        return pos
    return pos + max(1, svlen) - 1


def infer_svlen(ref, alt, info):
    svlen = as_int(info.get("SVLEN"))
    if svlen is not None:
        return svlen
    svtype = infer_svtype(ref, alt, info)
    if svtype == "DEL":
        return -abs(len(ref) - len(alt))
    if svtype == "INS":
        return abs(len(alt) - len(ref))
    end = as_int(info.get("END"))
    if end:
        return end
    return 0


def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def iter_vcf_records(path, chrom=None, start=None, end=None):
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            rec_chrom = parts[0]
            pos = int(parts[1])
            rec_id = parts[2]
            ref = parts[3]
            alt = parts[4]
            filt = parts[6]
            info = parse_info(parts[7])
            rec_svtype = infer_svtype(ref, alt, info)
            rec_end = infer_end(pos, ref, alt, info)
            rec_svlen = infer_svlen(ref, alt, info)
            if chrom and rec_chrom != chrom:
                continue
            if start is not None and rec_end < start:
                continue
            if end is not None and pos > end:
                continue
            yield {
                "chrom": rec_chrom,
                "pos": pos,
                "end": rec_end,
                "id": rec_id,
                "ref": ref,
                "alt": alt,
                "filter": filt,
                "info": info,
                "svtype": rec_svtype,
                "svlen": rec_svlen,
                "line": line.rstrip("\n"),
            }
