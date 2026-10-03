#!/usr/bin/env python3
import argparse
import bisect
import csv
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

def main():
    ap = argparse.ArgumentParser(description="Summarize DEL/INS panSV PAV overlap with Msa EDTA TE annotation using minimum overlap fraction.")
    ap.add_argument("--pav", required=True)
    ap.add_argument("--te-bed", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--exclude-sample", default="genome_Msa")
    ap.add_argument("--ins-window", type=int, default=100)
    ap.add_argument("--min-overlap-fraction", type=float, default=0.5)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    te_index = load_te_bed(args.te_bed)

    detailed = Counter(); plot3 = Counter(); totals = Counter()
    skipped_no_presence = 0; skipped_non_delins = 0
    overlap_fraction_bins = Counter()

    with open(args.pav, "r", encoding="utf-8", errors="replace", newline="") as inf:
        reader = csv.DictReader(inf, delimiter="\t")
        sample_cols = [c for c in reader.fieldnames if c not in META_COLS and c != args.exclude_sample]
        for row in reader:
            svtype = row.get("SVTYPE", "")
            if svtype not in ("DEL", "INS"):
                skipped_non_delins += 1
                continue
            present = sum(1 for c in sample_cols if presence_value(row.get(c)))
            if present == 0:
                skipped_no_presence += 1
                continue
            chrom, qstart, qend = sv_interval(row, args.ins_window)
            qlen = max(1, qend - qstart + 1)
            ov, union_ov = te_overlaps(te_index, chrom, qstart, qend)
            frac = min(1.0, union_ov / qlen)
            dom = dominant_class(ov)
            totals[svtype] += 1
            if union_ov == 0:
                detail = "Non-TE"
                plot_cls = "Non-TE"
                overlap_fraction_bins[(svtype, "0")] += 1
            elif frac < args.min_overlap_fraction:
                detail = "Low_overlap_TE"
                plot_cls = "Non-TE"
                overlap_fraction_bins[(svtype, "0-0.5") if args.min_overlap_fraction == 0.5 else (svtype, "below_threshold")] += 1
            else:
                detail = dom
                plot_cls = dom if dom in ("DNA", "LTR") else "Non-TE"
                overlap_fraction_bins[(svtype, ">=0.5") if args.min_overlap_fraction == 0.5 else (svtype, "pass_threshold")] += 1
            detailed[(svtype, detail)] += 1
            plot3[(svtype, plot_cls)] += 1

    detailed_out = args.out_prefix + ".detailed_counts.tsv"
    with open(detailed_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["SVTYPE", "TE_class_detailed", "count", "percent_of_svtype"])
        for svtype in ("DEL", "INS"):
            for cls in ("DNA", "LTR", "Other_TE", "Low_overlap_TE", "Non-TE"):
                cnt = detailed[(svtype, cls)]
                pct = cnt / totals[svtype] * 100 if totals[svtype] else 0
                w.writerow([svtype, cls, cnt, f"{pct:.4f}"])

    plot_out = args.out_prefix + ".plot3class.tsv"
    with open(plot_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["SVTYPE", "TE_class", "count", "percent_of_svtype"])
        for svtype in ("DEL", "INS"):
            for cls in ("DNA", "LTR", "Non-TE"):
                cnt = plot3[(svtype, cls)]
                pct = cnt / totals[svtype] * 100 if totals[svtype] else 0
                w.writerow([svtype, cls, cnt, f"{pct:.4f}"])

    frac_out = args.out_prefix + ".overlap_fraction_bins.tsv"
    with open(frac_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["SVTYPE", "overlap_fraction_bin", "count", "percent_of_svtype"])
        bins = ["0", "0-0.5", ">=0.5"] if args.min_overlap_fraction == 0.5 else ["0", "below_threshold", "pass_threshold"]
        for svtype in ("DEL", "INS"):
            for b in bins:
                cnt = overlap_fraction_bins[(svtype, b)]
                pct = cnt / totals[svtype] * 100 if totals[svtype] else 0
                w.writerow([svtype, b, cnt, f"{pct:.4f}"])

    summary_out = args.out_prefix + ".summary.tsv"
    with open(summary_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerow(["input_pav", args.pav])
        w.writerow(["te_bed", args.te_bed])
        w.writerow(["excluded_sample", args.exclude_sample])
        w.writerow(["remaining_sample_count", 17])
        w.writerow(["ins_breakpoint_window_bp", args.ins_window])
        w.writerow(["min_overlap_fraction", args.min_overlap_fraction])
        w.writerow(["DEL_total_used", totals["DEL"]])
        w.writerow(["INS_total_used", totals["INS"]])
        w.writerow(["skipped_non_DEL_INS", skipped_non_delins])
        w.writerow(["skipped_DEL_INS_no_remaining_presence", skipped_no_presence])
        w.writerow(["event_table", "not_written_to_avoid_large_intermediate_file"])
        w.writerow(["classification_note", "TE-associated requires union TE overlap/query interval >= min_overlap_fraction. DEL query interval is full Msa-reference SV interval; INS query interval is Msa-reference breakpoint +/- window proxy, not inserted-sequence TE annotation."])
    print(f"DONE DEL={totals['DEL']} INS={totals['INS']} skipped_no_presence={skipped_no_presence}")
    print(summary_out); print(plot_out); print(detailed_out); print(frac_out)

if __name__ == "__main__":
    main()
