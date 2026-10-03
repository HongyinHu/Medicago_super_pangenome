#!/usr/bin/env python3
import argparse
import csv
import os
from collections import Counter, defaultdict


MODELS = [
    "main_high_quality",
    "cladeI_internal_high_quality",
    "include_low_quality",
    "include_weak_spiny",
]


def read_candidates(path):
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = {}
        for row in reader:
            sv_id = row.get("MetaSV_ID") or row.get("SV_ID") or row.get("ID")
            if not sv_id:
                continue
            direction = row["association_direction"]
            rows[(sv_id, direction)] = row
    return rows


def load_dataset(dataset_dir):
    by_model = {}
    for model in MODELS:
        by_model[model] = read_candidates(os.path.join(dataset_dir, model, "candidates.tsv"))
    return by_model


def classify(models):
    main = "main_high_quality" in models
    clade = "cladeI_internal_high_quality" in models
    low = "include_low_quality" in models
    weak = "include_weak_spiny" in models
    if main and clade and low:
        return "Tier1_main_cladeI_low_quality_robust"
    if main and clade:
        return "Tier2_main_cladeI"
    if clade and low:
        return "Tier3_cladeI_low_quality_no_global"
    if main:
        return "Tier4_main_only_phylo_sensitive"
    if clade:
        return "Tier5_cladeI_only"
    if low or weak:
        return "Tier6_sensitivity_only"
    return "unclassified"


def tier_dataset(dataset, dataset_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    by_model = load_dataset(dataset_dir)
    keys = set()
    for rows in by_model.values():
        keys.update(rows)

    out_rows = []
    for key in sorted(keys):
        models = [model for model in MODELS if key in by_model[model]]
        representative = None
        for model in MODELS:
            if key in by_model[model]:
                representative = dict(by_model[model][key])
                break
        representative["support_models"] = ",".join(models)
        representative["n_support_models"] = str(len(models))
        representative["tier"] = classify(set(models))
        representative["main_support"] = "1" if "main_high_quality" in models else "0"
        representative["cladeI_internal_support"] = "1" if "cladeI_internal_high_quality" in models else "0"
        representative["include_low_quality_support"] = "1" if "include_low_quality" in models else "0"
        representative["include_weak_spiny_support"] = "1" if "include_weak_spiny" in models else "0"
        out_rows.append(representative)

    field_order = [
        "dataset",
        "tier",
        "support_models",
        "n_support_models",
        "main_support",
        "cladeI_internal_support",
        "include_low_quality_support",
        "include_weak_spiny_support",
        "association_direction",
    ]
    if out_rows:
        for key in out_rows[0]:
            if key not in field_order:
                field_order.append(key)

    tiered_path = os.path.join(out_dir, "tiered_candidates.tsv")
    with open(tiered_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=field_order)
        writer.writeheader()
        writer.writerows(out_rows)

    high_conf = [
        row
        for row in out_rows
        if row["tier"] in {"Tier1_main_cladeI_low_quality_robust", "Tier2_main_cladeI"}
    ]
    high_conf_path = os.path.join(out_dir, "high_confidence_candidates.tsv")
    with open(high_conf_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=field_order)
        writer.writeheader()
        writer.writerows(high_conf)

    summary = Counter((row["tier"], row["association_direction"]) for row in out_rows)
    svtype_summary = Counter((row["tier"], row.get("SVTYPE", ".")) for row in out_rows)
    with open(os.path.join(out_dir, "tier_summary.tsv"), "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["dataset", "section", "tier", "key", "count"])
        for (tier, direction), count in sorted(summary.items()):
            writer.writerow([dataset, "direction", tier, direction, count])
        for (tier, svtype), count in sorted(svtype_summary.items()):
            writer.writerow([dataset, "svtype", tier, svtype, count])
    return out_rows, high_conf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True)
    args = parser.parse_args()

    datasets = ["meta_publication", "RefA", "RefB"]
    combined = []
    combined_high = []
    for dataset in datasets:
        rows, high = tier_dataset(
            dataset,
            os.path.join(args.results_dir, dataset),
            os.path.join(args.results_dir, dataset),
        )
        combined.extend(rows)
        combined_high.extend(high)

    if combined:
        fieldnames = []
        preferred = [
            "dataset",
            "tier",
            "support_models",
            "n_support_models",
            "main_support",
            "cladeI_internal_support",
            "include_low_quality_support",
            "include_weak_spiny_support",
            "association_direction",
            "MetaSV_ID",
            "SV_ID",
            "RefA_SV_ID",
            "RefB_SV_ID",
            "CHROM",
            "POS",
            "END",
            "SVTYPE",
            "SVLEN",
        ]
        all_keys = set()
        for row in combined:
            all_keys.update(row)
        for key in preferred:
            if key in all_keys:
                fieldnames.append(key)
        for key in sorted(all_keys):
            if key not in fieldnames:
                fieldnames.append(key)
        with open(os.path.join(args.results_dir, "all_datasets_tiered_candidates.tsv"), "w", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(combined)
        with open(os.path.join(args.results_dir, "all_datasets_high_confidence_candidates.tsv"), "w", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(combined_high)


if __name__ == "__main__":
    main()
