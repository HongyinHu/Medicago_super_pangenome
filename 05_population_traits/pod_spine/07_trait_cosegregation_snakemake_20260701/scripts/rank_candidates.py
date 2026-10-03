#!/usr/bin/env python3
import argparse
import csv
import os


def read_optional(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--event-level", required=True)
    parser.add_argument("--gene-interval", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    rows = []
    for path in (args.event_level, args.gene_interval):
        for row in read_optional(path):
            pos = fnum(row.get("positive_match_rate"))
            neg = fnum(row.get("negative_match_rate"))
            pass_flag = int(row.get("pass_80pct", "0") or 0)
            row["group_balance_score"] = f"{min(pos, neg):.3f}"
            row["priority_score"] = f"{(100 * pass_flag + 10 * min(pos, neg)):.3f}"
            rows.append(row)

    rows.sort(key=lambda r: fnum(r.get("priority_score")), reverse=True)
    field_order = [
        "candidate_id", "evidence_track", "priority_score", "group_balance_score",
        "pass_80pct", "priority_note", "ref", "chrom", "start", "end", "gene_id",
        "region_type", "expected_state", "best_direction",
        "positive_match", "positive_total", "positive_match_rate",
        "negative_match", "negative_total", "negative_match_rate",
    ]
    fields = field_order + sorted({k for r in rows for k in r if k not in field_order})
    with open(args.out, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
