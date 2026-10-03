#!/usr/bin/env python3
import argparse
import csv
import os


PRESENT = {"1", "present", "P", "TRUE", "True", "true", "0/1", "1/1"}
ABSENT = {"0", "absent", "A", "FALSE", "False", "false", "0/0"}
MISSING = {".", "./.", "NA", "na", "NaN", "", "missing"}


def read_phenotype(path):
    pheno = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row["include_primary"] = str(row.get("include_primary", "1")) == "1"
            pheno[row["sample"]] = row
    return pheno


def normalize_state(value):
    value = str(value).strip()
    if value in PRESENT:
        return "present"
    if value in ABSENT:
        return "absent"
    if value in MISSING:
        return "missing"
    return "missing"


def score_event(row, samples, pheno, positive, negative, expected_present_in):
    pos_total = pos_match = neg_total = neg_match = missing = 0
    for sample in samples:
        meta = pheno.get(sample)
        if not meta or not meta.get("include_primary"):
            continue
        state = normalize_state(row.get(sample, "missing"))
        if state == "missing":
            missing += 1
            continue
        phenotype = meta["phenotype"]
        expect_present = phenotype == expected_present_in
        match = (state == "present") if expect_present else (state == "absent")
        if phenotype == positive:
            pos_total += 1
            pos_match += int(match)
        elif phenotype == negative:
            neg_total += 1
            neg_match += int(match)
    pos_rate = pos_match / pos_total if pos_total else 0.0
    neg_rate = neg_match / neg_total if neg_total else 0.0
    return pos_total, pos_match, pos_rate, neg_total, neg_match, neg_rate, missing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--pav", required=True)
    parser.add_argument("--positive", required=True)
    parser.add_argument("--negative", required=True)
    parser.add_argument("--min-match", type=float, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pheno = read_phenotype(args.phenotype)

    with open(args.pav, newline="", encoding="utf-8") as handle, open(args.out, "w", newline="", encoding="utf-8") as out:
        reader = csv.DictReader(handle, delimiter="\t")
        samples = [c for c in reader.fieldnames if c != "event_id"]
        fields = [
            "candidate_id", "evidence_track", "best_direction",
            "positive_match", "positive_total", "positive_match_rate",
            "negative_match", "negative_total", "negative_match_rate",
            "missing_primary", "pass_80pct", "priority_note",
        ]
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in reader:
            scores = []
            for direction in (args.positive, args.negative):
                scores.append((direction,) + score_event(row, samples, pheno, args.positive, args.negative, direction))
            best = max(scores, key=lambda x: min(x[3], x[6]))
            direction, pos_total, pos_match, pos_rate, neg_total, neg_match, neg_rate, missing = best
            pass_flag = pos_rate >= args.min_match and neg_rate >= args.min_match
            writer.writerow({
                "candidate_id": row.get("event_id", ""),
                "evidence_track": "event_level_pav",
                "best_direction": f"present_in_{direction}",
                "positive_match": pos_match,
                "positive_total": pos_total,
                "positive_match_rate": f"{pos_rate:.3f}",
                "negative_match": neg_match,
                "negative_total": neg_total,
                "negative_match_rate": f"{neg_rate:.3f}",
                "missing_primary": missing,
                "pass_80pct": int(pass_flag),
                "priority_note": "event_level_pass" if pass_flag else "event_level_fail",
            })


if __name__ == "__main__":
    main()
