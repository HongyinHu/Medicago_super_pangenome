#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from bisect import bisect_right
from collections import Counter, defaultdict


SKIP_FEATURES = {"target_site_duplication"}


def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_attrs(attr):
    out = {}
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.lower()] = v
        elif " " in part:
            k, v = part.split(" ", 1)
            out[k.lower()] = v.strip('"')
    return out


def load_te(gff):
    te = defaultdict(list)
    with open_maybe_gzip(gff) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = parts
            if feature in SKIP_FEATURES:
                continue
            a = parse_attrs(attrs)
            classification = a.get("classification") or a.get("class") or feature
            name = a.get("name") or a.get("id") or feature
            te[chrom].append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "feature": feature,
                    "class": classification,
                    "name": name,
                }
            )
    index = {}
    for chrom, arr in te.items():
        arr.sort(key=lambda x: (x["start"], x["end"]))
        starts = [x["start"] for x in arr]
        index[chrom] = (arr, starts)
    return index


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def annotate_interval(index, chrom, start, end):
    if not chrom or chrom in {"NA", "."} or not start or not end:
        return {"bp": 0, "fraction": 0.0, "classes": "NA", "features": "NA", "names": "NA", "count": 0}
    try:
        start, end = int(start), int(end)
    except ValueError:
        return {"bp": 0, "fraction": 0.0, "classes": "NA", "features": "NA", "names": "NA", "count": 0}
    if start > end:
        start, end = end, start
    length = max(1, end - start + 1)
    if chrom not in index:
        return {"bp": 0, "fraction": 0.0, "classes": "chrom_missing_in_te", "features": "NA", "names": "NA", "count": 0}
    arr, starts = index[chrom]
    pos = bisect_right(starts, end)
    scan_start = max(0, pos - 1000)
    overlaps = []
    classes = Counter()
    features = Counter()
    names = Counter()
    for te in arr[scan_start:pos]:
        if te["end"] < start:
            continue
        if te["start"] > end:
            break
        s = max(start, te["start"])
        e = min(end, te["end"])
        if s <= e:
            overlaps.append((s, e))
            classes[te["class"]] += 1
            features[te["feature"]] += 1
            names[te["name"]] += 1
    merged = merge_intervals(overlaps)
    bp = sum(e - s + 1 for s, e in merged)
    frac = min(1.0, bp / length)
    return {
        "bp": bp,
        "fraction": round(frac, 4),
        "classes": ",".join(k for k, _ in classes.most_common(8)) if classes else "NA",
        "features": ",".join(k for k, _ in features.most_common(8)) if features else "NA",
        "names": ",".join(k for k, _ in names.most_common(8)) if names else "NA",
        "count": sum(classes.values()),
    }


def as_float(v):
    try:
        return float(v)
    except Exception:
        return 0.0


def te_flag(frac):
    if frac >= 0.9:
        return "fully_TE_embedded"
    if frac >= 0.5:
        return "high_TE_overlap"
    if frac >= 0.1:
        return "partial_TE_overlap"
    if frac > 0:
        return "minor_TE_overlap"
    return "no_TE_overlap"


def te_score_adjust(flag):
    if flag == "no_TE_overlap":
        return 10
    if flag == "minor_TE_overlap":
        return 5
    if flag == "partial_TE_overlap":
        return 0
    if flag == "high_TE_overlap":
        return -6
    if flag == "fully_TE_embedded":
        return -10
    return 0


def annotate(input_path, output_path, refa_te, refb_te):
    indexes = {"RefA": load_te(refa_te), "RefB": load_te(refb_te)}
    rows = []
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        base_fields = reader.fieldnames or []
        for row in reader:
            best_frac = 0.0
            best_ref = "NA"
            best_classes = "NA"
            for ref in ["RefA", "RefB"]:
                anno = annotate_interval(indexes[ref], row.get(f"{ref}_CHROM", ""), row.get(f"{ref}_POS", ""), row.get(f"{ref}_END", ""))
                row[f"{ref}_TE_overlap_bp"] = anno["bp"]
                row[f"{ref}_TE_overlap_fraction"] = anno["fraction"]
                row[f"{ref}_TE_overlap_count"] = anno["count"]
                row[f"{ref}_TE_classes"] = anno["classes"]
                row[f"{ref}_TE_features"] = anno["features"]
                row[f"{ref}_TE_names"] = anno["names"]
                frac = as_float(anno["fraction"])
                if frac > best_frac:
                    best_frac = frac
                    best_ref = ref
                    best_classes = anno["classes"]
            flag = te_flag(best_frac)
            row["best_TE_ref"] = best_ref
            row["best_TE_overlap_fraction"] = round(best_frac, 4)
            row["best_TE_classes"] = best_classes
            row["TE_overlap_flag"] = flag
            row["TE_score_adjustment"] = te_score_adjust(flag)
            base_score = as_float(row.get("gene_validation_score") or row.get("validation_score"))
            row["final_validation_score"] = round(base_score + te_score_adjust(flag), 3)
            rows.append(row)

    rows.sort(
        key=lambda r: (
            -as_float(r["final_validation_score"]),
            -as_float(r.get("gene_validation_score") or r.get("validation_score")),
            r.get("TE_overlap_flag", ""),
            r.get("MetaSV_ID", ""),
        )
    )
    for i, row in enumerate(rows, 1):
        row["final_validation_rank"] = i

    added = [
        "final_validation_rank",
        "final_validation_score",
        "TE_score_adjustment",
        "TE_overlap_flag",
        "best_TE_ref",
        "best_TE_overlap_fraction",
        "best_TE_classes",
        "RefA_TE_overlap_bp",
        "RefA_TE_overlap_fraction",
        "RefA_TE_overlap_count",
        "RefA_TE_classes",
        "RefA_TE_features",
        "RefA_TE_names",
        "RefB_TE_overlap_bp",
        "RefB_TE_overlap_fraction",
        "RefB_TE_overlap_count",
        "RefB_TE_classes",
        "RefB_TE_features",
        "RefB_TE_names",
    ]
    fields = added + [f for f in base_fields if f not in added]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_subset(path, rows, predicate):
    with open(path, "w", newline="") as f:
        fields = rows[0].keys() if rows else []
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            if predicate(row):
                writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--refa-te", required=True)
    parser.add_argument("--refb-te", required=True)
    args = parser.parse_args()
    rows = annotate(args.input, args.output, args.refa_te, args.refb_te)
    prefix = re.sub(r"\.tsv$", "", args.output)
    write_subset(prefix + ".low_TE_first_review.tsv", rows, lambda r: r["TE_overlap_flag"] in {"no_TE_overlap", "minor_TE_overlap"} and r.get("priority_class") == "A_top_manual_review")
    write_subset(prefix + ".TE_related.tsv", rows, lambda r: r["TE_overlap_flag"] in {"high_TE_overlap", "fully_TE_embedded"} and r.get("priority_class") in {"A_top_manual_review", "B_strong_candidate"})
    print(f"wrote\t{args.output}\trows={len(rows)}")


if __name__ == "__main__":
    main()
