#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def as_int(value, default=0):
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except ValueError:
        return default


def as_float(value, default=0.0):
    try:
        if value is None or value == "":
            return default
        return float(value)
    except ValueError:
        return default


def bool_text(flag):
    return "TRUE" if flag else "FALSE"


def length_score(abs_len):
    if 50 <= abs_len <= 5_000:
        return 12, "PCR_friendly_50bp_5kb"
    if 5_000 < abs_len <= 20_000:
        return 7, "moderate_5kb_20kb"
    if 20_000 < abs_len <= 100_000:
        return 2, "large_20kb_100kb"
    return -8, "very_large_gt100kb"


def classify_priority(score, exact, tier, cross_ref):
    if score >= 120 or (exact and tier.startswith("Tier1") and score >= 105):
        return "A_top_manual_review"
    if score >= 95:
        return "B_strong_candidate"
    if score >= 75:
        return "C_secondary_candidate"
    return "D_low_priority_within_high_conf"


def score_row(row):
    tier = row.get("tier", "")
    direction = row.get("association_direction", "")
    svtype = row.get("SVTYPE", "")
    match_status = row.get("match_status", "")
    confidence = row.get("confidence", "")

    spiny_called = as_int(row.get("spiny_called"))
    spiny_present = as_int(row.get("spiny_present"))
    spiny_absent = as_int(row.get("spiny_absent"))
    spiny_missing = as_int(row.get("spiny_missing"))
    spineless_called = as_int(row.get("spineless_called"))
    spineless_present = as_int(row.get("spineless_present"))
    spineless_absent = as_int(row.get("spineless_absent"))
    spineless_missing = as_int(row.get("spineless_missing"))
    spiny_rate = as_float(row.get("spiny_present_rate"))
    spineless_rate = as_float(row.get("spineless_present_rate"))

    score = 0.0
    reasons = []

    if tier.startswith("Tier1"):
        score += 35
        reasons.append("Tier1 robust to low-quality sensitivity")
    elif tier.startswith("Tier2"):
        score += 22
        reasons.append("Tier2 main plus Clade-I support")

    if row.get("cladeI_internal_support", "").upper() == "TRUE":
        score += 15
        reasons.append("Clade-I internal support")
    if row.get("include_low_quality_support", "").upper() == "TRUE":
        score += 15
        reasons.append("robust after including low-quality genomes")
    if row.get("include_weak_spiny_support", "").upper() == "TRUE":
        score += 8
        reasons.append("robust when weak-spiny genome_436 is included")

    if direction == "spiny_present":
        phenotype_error_count = spiny_absent + spineless_present
        exact = (spiny_absent == 0 and spineless_present == 0)
        separation = max(0.0, spiny_rate - spineless_rate)
    else:
        phenotype_error_count = spiny_present + spineless_absent
        exact = (spiny_present == 0 and spineless_absent == 0)
        separation = max(0.0, spineless_rate - spiny_rate)

    score += separation * 28
    if exact:
        score += 28
        reasons.append("no phenotype-discordant called samples")
    else:
        score -= phenotype_error_count * 12
        reasons.append(f"{phenotype_error_count} discordant called samples")

    missing_total = spiny_missing + spineless_missing
    score -= missing_total * 2
    if spiny_called >= 4:
        score += 4
    if spineless_called >= 9:
        score += 4

    svtype_score = {"DEL": 12, "INS": 10, "INV": 3, "DUP": 3}.get(svtype, -4)
    score += svtype_score
    if svtype in {"DEL", "INS"}:
        reasons.append(f"{svtype} is validation-friendly")
    else:
        reasons.append(f"{svtype} is lower priority for PCR validation")

    abs_len = abs(as_int(row.get("SVLEN")))
    l_score, l_class = length_score(abs_len)
    score += l_score
    reasons.append(l_class)

    cross_ref_supported = match_status.startswith("shared_strict")
    if cross_ref_supported:
        score += 20
        reasons.append("cross-reference strict support")
    else:
        reasons.append("reference-specific candidate")

    if confidence == "high":
        score += 8
        reasons.append("high confidence")

    priority_class = classify_priority(score, exact, tier, cross_ref_supported)
    if exact and svtype in {"DEL", "INS"} and abs_len <= 20_000:
        validation_bucket = "first_batch"
    elif svtype in {"DEL", "INS"} and phenotype_error_count <= 1 and abs_len <= 100_000:
        validation_bucket = "second_batch"
    else:
        validation_bucket = "deprioritize_or_visual_only"

    return {
        "validation_score": round(score, 3),
        "priority_class": priority_class,
        "validation_bucket": validation_bucket,
        "exact_cosegregation": bool_text(exact),
        "phenotype_error_count": phenotype_error_count,
        "missing_total": missing_total,
        "abs_svlen": abs_len,
        "size_class": l_class,
        "cross_ref_supported": bool_text(cross_ref_supported),
        "separation_delta": round(separation, 4),
        "key_reason": "; ".join(reasons),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--top", type=int, default=50)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = []
        for row in reader:
            metrics = score_row(row)
            rows.append({**metrics, **row})
        original_fields = reader.fieldnames or []

    rows.sort(
        key=lambda r: (
            -as_float(r["validation_score"]),
            r["priority_class"],
            r.get("tier", ""),
            r.get("SVTYPE", ""),
            as_int(r.get("abs_svlen")),
            r.get("MetaSV_ID", ""),
        )
    )
    for i, row in enumerate(rows, 1):
        row["validation_rank"] = i

    new_fields = [
        "validation_rank",
        "validation_score",
        "priority_class",
        "validation_bucket",
        "exact_cosegregation",
        "phenotype_error_count",
        "missing_total",
        "abs_svlen",
        "size_class",
        "cross_ref_supported",
        "separation_delta",
        "key_reason",
    ]
    fields = new_fields + [f for f in original_fields if f not in new_fields]

    def write_table(path, subset):
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(subset)

    write_table(outdir / "validation_priority_candidates.tsv", rows)
    write_table(outdir / f"validation_priority_top{args.top}.tsv", rows[: args.top])
    write_table(outdir / "validation_priority_top100.tsv", rows[:100])
    write_table(outdir / "validation_priority_first_batch.tsv", [r for r in rows if r["validation_bucket"] == "first_batch"])
    write_table(outdir / "validation_priority_A_candidates.tsv", [r for r in rows if r["priority_class"] == "A_top_manual_review"])

    counters = {}
    for row in rows:
        for name, key in [
            ("priority_class", row["priority_class"]),
            ("validation_bucket", row["validation_bucket"]),
            ("tier", row.get("tier", "")),
            ("direction", row.get("association_direction", "")),
            ("svtype", row.get("SVTYPE", "")),
            ("cross_ref_supported", row["cross_ref_supported"]),
            ("exact_cosegregation", row["exact_cosegregation"]),
        ]:
            counters[(name, key)] = counters.get((name, key), 0) + 1

    with open(outdir / "validation_priority_summary.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["category", "value", "count"])
        writer.writerow(["total", "all", len(rows)])
        for (category, value), count in sorted(counters.items()):
            writer.writerow([category, value, count])


if __name__ == "__main__":
    main()
