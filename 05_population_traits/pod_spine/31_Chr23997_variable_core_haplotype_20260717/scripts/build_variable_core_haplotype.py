#!/usr/bin/env python3
"""Summarize a variable-length Chr23997 intron-2 core-loss haplotype."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def spans_shared_core(pos, end, core_start, core_end, tolerance=20, min_length=80):
    """Return whether a DEL of any length spans the pre-defined shared core."""
    return (
        int(end) - int(pos) >= min_length
        and int(pos) <= core_start + tolerance
        and int(end) >= core_end - tolerance
    )


def is_local_core_deletion(row, core_start, core_end, tolerance=20, min_length=80, max_length=2000):
    """Reject chromosome-scale raw calls before scoring a local intronic event."""
    if row.get("svtype") != "DEL":
        return False
    length = abs(int(float(row.get("svlen") or 0)))
    return (
        min_length <= length <= max_length
        and spans_shared_core(row["pos"], row["end"], core_start, core_end, tolerance, min_length)
    )


def truth_state(value):
    """Keep the prior IGV-calibrated truth call separate from raw caller noise."""
    if value == "DEL_present":
        return "CORE_LOSS"
    if value == "INTACT_absent":
        return "INTACT"
    return "EXCLUDE"


def phenotype_from_group(group):
    group = group.lower()
    if "spineless" in group:
        return "spineless"
    if "spiny" in group:
        return "spiny"
    return ""


def build_primary_panel(rows, excluded_samples):
    """Select the read-validated species panel without forcing exceptions away."""
    panel = []
    for row in rows:
        if row["sample"] in excluded_samples:
            continue
        state = truth_state(row["validated_truth_class"])
        phenotype = phenotype_from_group(row["phenotype_group_curated"])
        if state == "EXCLUDE" or not phenotype:
            continue
        panel.append({"sample": row["sample"], "phenotype": phenotype, "core_state": state})
    return panel


def fisher_two_sided(a, b, c, d):
    """Two-sided Fisher exact test using the probability ordering definition."""
    row1, row2 = a + b, c + d
    col1, total = a + c, a + b + c + d

    def prob(x):
        return math.comb(col1, x) * math.comb(total - col1, row1 - x) / math.comb(total, row1)

    observed = prob(a)
    low = max(0, row1 - (total - col1))
    high = min(row1, col1)
    return sum(prob(x) for x in range(low, high + 1) if prob(x) <= observed + 1e-12)


def odds_ratio(a, b, c, d):
    denominator = b * c
    if denominator == 0:
        return "inf" if a * d else "NA"
    return f"{a * d / denominator:.6g}"


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-table", required=True)
    parser.add_argument("--caller-records", required=True)
    parser.add_argument("--out-haplotypes", required=True)
    parser.add_argument("--out-association", required=True)
    parser.add_argument("--core-start", type=int, default=88980879)
    parser.add_argument("--core-end", type=int, default=88981002)
    parser.add_argument("--tolerance", type=int, default=20)
    parser.add_argument("--max-local-sv-length", type=int, default=2000)
    parser.add_argument("--exclude-samples", default="genome_A17,genome_M46,genome_482,genome_Msa")
    args = parser.parse_args()

    excluded = {sample for sample in args.exclude_samples.split(",") if sample}
    raw_support = defaultdict(list)
    for row in read_tsv(args.caller_records):
        if row.get("filter") not in {"PASS", "."}:
            continue
        try:
            spans = is_local_core_deletion(
                row, args.core_start, args.core_end, args.tolerance, max_length=args.max_local_sv_length
            )
        except (KeyError, ValueError):
            continue
        if spans:
            raw_support[row["sample"]].append(row)

    truth_rows = read_tsv(args.truth_table)
    output = []
    for row in truth_rows:
        sample = row["sample"]
        calls = raw_support.get(sample, [])
        callers = sorted({call["caller"] for call in calls})
        events = ";".join(
            f"{call['caller']}:{call['pos']}-{call['end']}({call.get('svlen', '.')})" for call in calls
        ) or "NA"
        state = truth_state(row["validated_truth_class"])
        phenotype = phenotype_from_group(row["phenotype_group_curated"])
        output.append({
            "sample": sample,
            "phenotype_group_curated": row["phenotype_group_curated"],
            "phenotype": phenotype or "NA",
            "validated_truth_class": row["validated_truth_class"],
            "core_state": state,
            "in_primary_panel": "yes" if sample not in excluded and state != "EXCLUDE" and phenotype else "no",
            "raw_core_callers": ",".join(callers) or "NA",
            "n_raw_core_callers": len(callers),
            "raw_variable_length_events": events,
            "manual_igv_note": row.get("note", ""),
        })

    fields = [
        "sample", "phenotype_group_curated", "phenotype", "validated_truth_class", "core_state",
        "in_primary_panel", "raw_core_callers", "n_raw_core_callers", "raw_variable_length_events",
        "manual_igv_note",
    ]
    write_tsv(args.out_haplotypes, output, fields)

    panel = build_primary_panel(truth_rows, excluded)
    counts = defaultdict(int)
    for row in panel:
        counts[(row["phenotype"], row["core_state"])] += 1
    a = counts[("spineless", "CORE_LOSS")]
    b = counts[("spineless", "INTACT")]
    c = counts[("spiny", "CORE_LOSS")]
    d = counts[("spiny", "INTACT")]
    summary = [{
        "analysis": "read_validated_species_panel_unadjusted",
        "n_species": len(panel),
        "spineless_core_loss": a,
        "spineless_intact": b,
        "spiny_core_loss": c,
        "spiny_intact": d,
        "odds_ratio_spineless_core_loss": odds_ratio(a, b, c, d),
        "fisher_two_sided_p": f"{fisher_two_sided(a, b, c, d):.8g}",
        "scope": "Descriptive species-panel Fisher test; phylogenetic dependence is not adjusted.",
    }]
    write_tsv(args.out_association, summary, list(summary[0]))


if __name__ == "__main__":
    main()
