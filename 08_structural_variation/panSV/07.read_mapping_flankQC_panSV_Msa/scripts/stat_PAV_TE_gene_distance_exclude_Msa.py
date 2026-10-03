#!/usr/bin/env python3
import argparse
import bisect
import csv
import math
import os
from collections import defaultdict, Counter

META_COLS = {
    "SV_ID", "source_methods", "reference", "CHROM", "POS", "END", "SVTYPE", "SVLEN",
    "read_SV_ID", "svgap_ids", "final_confidence"
}

def presence_value(v):
    if v is None:
        return False
    s = str(v).strip()
    if s in ("", ".", "NA", "NaN", "nan"):
        return False
    try:
        return float(s) > 0
    except ValueError:
        return s not in ("0", "./.", "0/0")

def classify_te(fields):
    text = "\t".join(fields).lower()
    if "ltr" in text:
        return "LTR"
    if "helitron" in text or "tir" in text or "mite" in text or "/dt" in text or "dna" in text:
        return "DNA"
    return "Other_TE"

def load_te_bed(te_bed):
    intervals = defaultdict(list)
    with open(te_bed, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start = int(float(f[1])); end = int(float(f[2]))
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            intervals[f[0]].append((start, end, classify_te(f)))
    index = {}
    for chrom, recs in intervals.items():
        recs.sort(key=lambda x: (x[0], x[1]))
        starts = [r[0] for r in recs]
        prefix = []
        m = -1
        for r in recs:
            m = max(m, r[1]); prefix.append(m)
        index[chrom] = {"recs": recs, "starts": starts, "prefix": prefix}
    return index

def merge_len(intervals):
    if not intervals:
        return 0
    intervals.sort()
    total = 0
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s, e
    total += ce - cs + 1
    return total

def te_overlaps(te_index, chrom, start, end):
    idx = te_index.get(chrom)
    if not idx:
        return Counter(), 0
    recs, starts, prefix = idx["recs"], idx["starts"], idx["prefix"]
    j = bisect.bisect_right(starts, end) - 1
    counts = Counter(); union = []
    while j >= 0 and prefix[j] >= start:
        ts, te, cls = recs[j]
        if te >= start and ts <= end:
            os, oe = max(start, ts), min(end, te)
            if oe >= os:
                counts[cls] += oe - os + 1
                union.append((os, oe))
        j -= 1
    return counts, merge_len(union)

def dominant_class(counts):
    if not counts:
        return "Non-TE"
    order = {"LTR": 0, "DNA": 1, "Other_TE": 2}
    return sorted(counts.items(), key=lambda x: (-x[1], order.get(x[0], 9)))[0][0]

def sv_interval(row, ins_window):
    chrom = row["CHROM"]
    pos = int(float(row["POS"]))
    try:
        end = int(float(row.get("END", pos))) if row.get("END", "") not in ("", ".", "NA") else pos
    except ValueError:
        end = pos
    if end < pos:
        pos, end = end, pos
    if row.get("SVTYPE", "") == "INS":
        return chrom, max(1, pos - ins_window), pos + ins_window
    return chrom, pos, end

def quantile(vals, q):
    if not vals:
        return "."
    vals = sorted(vals)
    pos = (len(vals) - 1) * q
    lo = int(math.floor(pos)); hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)

def load_distance(gene_context_events):
    distance_by_id = {}
    with open(gene_context_events, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sv_id = row.get("SV_ID")
            if not sv_id:
                continue
            d = row.get("distance_to_gene_bp", ".")
            if d in ("", ".", "NA", "nan"):
                continue
            try:
                distance_by_id[sv_id] = int(float(d))
            except ValueError:
                continue
    return distance_by_id

def main():
    ap = argparse.ArgumentParser(description="Compute nearest-gene distance by strict TE-associated PAV class.")
    ap.add_argument("--pav", required=True)
    ap.add_argument("--te-bed", required=True)
    ap.add_argument("--gene-context-events", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--exclude-sample", default="genome_Msa")
    ap.add_argument("--ins-window", type=int, default=100)
    ap.add_argument("--min-overlap-fraction", type=float, default=0.5)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    te_index = load_te_bed(args.te_bed)
    distance_by_id = load_distance(args.gene_context_events)

    values_out = args.out_prefix + ".values.tsv"
    summary = defaultdict(list)
    detailed = Counter()
    totals = Counter()
    missing_distance = 0; skipped_no_presence = 0; skipped_non_delins = 0; used = 0

    with open(args.pav, "r", encoding="utf-8", errors="replace", newline="") as inf, \
         open(values_out, "w", encoding="utf-8", newline="") as outf:
        reader = csv.DictReader(inf, delimiter="\t")
        sample_cols = [c for c in reader.fieldnames if c not in META_COLS and c != args.exclude_sample]
        w = csv.writer(outf, delimiter="\t", lineterminator="\n")
        w.writerow(["TE_class_plot", "distance_to_gene_bp", "log2_distance_plus1"])
        for row in reader:
            svtype = row.get("SVTYPE", "")
            if svtype not in ("DEL", "INS"):
                skipped_non_delins += 1
                continue
            present = sum(1 for c in sample_cols if presence_value(row.get(c)))
            if present == 0:
                skipped_no_presence += 1
                continue
            sv_id = row.get("SV_ID")
            dist = distance_by_id.get(sv_id)
            if dist is None:
                missing_distance += 1
                continue
            chrom, qstart, qend = sv_interval(row, args.ins_window)
            qlen = max(1, qend - qstart + 1)
            ov, union_ov = te_overlaps(te_index, chrom, qstart, qend)
            frac = min(1.0, union_ov / qlen)
            dom = dominant_class(ov)
            if union_ov > 0 and frac >= args.min_overlap_fraction and dom in ("DNA", "LTR"):
                cls = dom
            else:
                cls = "Non-TE"
            totals[svtype] += 1
            detailed[(svtype, cls)] += 1
            logd = math.log2(dist + 1)
            w.writerow([cls, dist, f"{logd:.8f}"])
            summary[cls].append(dist)
            used += 1

    summary_out = args.out_prefix + ".summary.tsv"
    with open(summary_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["TE_class", "count", "median_distance_bp", "mean_distance_bp", "q1_distance_bp", "q3_distance_bp", "median_log2_distance_plus1", "mean_log2_distance_plus1"])
        for cls in ("Non-TE", "DNA", "LTR"):
            vals = summary.get(cls, [])
            logs = [math.log2(v + 1) for v in vals]
            if vals:
                w.writerow([cls, len(vals), f"{quantile(vals, 0.5):.2f}", f"{sum(vals)/len(vals):.2f}", f"{quantile(vals, 0.25):.2f}", f"{quantile(vals, 0.75):.2f}", f"{quantile(logs, 0.5):.4f}", f"{sum(logs)/len(logs):.4f}"])
            else:
                w.writerow([cls, 0, ".", ".", ".", ".", ".", "."])

    class_counts_out = args.out_prefix + ".class_counts.tsv"
    with open(class_counts_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["SVTYPE", "TE_class_plot", "count"])
        for svtype in ("DEL", "INS"):
            for cls in ("Non-TE", "DNA", "LTR"):
                w.writerow([svtype, cls, detailed[(svtype, cls)]])

    run_out = args.out_prefix + ".run_summary.tsv"
    with open(run_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerow(["pav", args.pav])
        w.writerow(["te_bed", args.te_bed])
        w.writerow(["gene_context_events", args.gene_context_events])
        w.writerow(["values_table", values_out])
        w.writerow(["excluded_sample", args.exclude_sample])
        w.writerow(["ins_breakpoint_window_bp", args.ins_window])
        w.writerow(["min_overlap_fraction", args.min_overlap_fraction])
        w.writerow(["events_used", used])
        w.writerow(["missing_distance", missing_distance])
        w.writerow(["skipped_non_DEL_INS", skipped_non_delins])
        w.writerow(["skipped_DEL_INS_no_remaining_presence", skipped_no_presence])
        w.writerow(["classification_note", "TE-associated classes require union TE overlap/query interval >= min_overlap_fraction; otherwise included in Non-TE for this distance comparison."])
    print(f"DONE events_used={used} missing_distance={missing_distance}")
    print(values_out); print(summary_out)

if __name__ == "__main__":
    main()
