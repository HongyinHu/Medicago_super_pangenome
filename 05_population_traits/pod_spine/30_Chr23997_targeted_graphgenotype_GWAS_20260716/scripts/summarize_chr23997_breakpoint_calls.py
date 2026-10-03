#!/usr/bin/env python3
"""Merge Chr23997 core-loss genotypes and report the unadjusted contingency test."""

import argparse
import csv
from pathlib import Path

from scipy.stats import fisher_exact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--calls-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    with open(args.manifest, newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for sample in samples:
        path = Path(args.calls_dir) / f"{sample['sample']}.tsv"
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            call = next(csv.DictReader(handle, delimiter="\t"))
        call.update({"phenotype": sample["phenotype"], "task_index": sample["task_index"]})
        call["del_present"] = "1" if call["GT"] in {"0/1", "1/1"} else "0"
        rows.append(call)

    fields = ["task_index", "sample", "phenotype", "GT", "del_present", "del_split_templates",
              "ref_continuous_templates", "flank_mean_depth", "core_mean_depth", "core_to_flank_ratio",
              "modal_del_breakpoint_1based", "evidence_class"]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda item: int(item["task_index"])))

    called = [row for row in rows if row["GT"] != "./."]
    case_del = sum(row["phenotype"] == "1" and row["del_present"] == "1" for row in called)
    case_ref = sum(row["phenotype"] == "1" and row["del_present"] == "0" for row in called)
    control_del = sum(row["phenotype"] == "0" and row["del_present"] == "1" for row in called)
    control_ref = sum(row["phenotype"] == "0" and row["del_present"] == "0" for row in called)
    odds_ratio, p_value = fisher_exact([[case_del, case_ref], [control_del, control_ref]], alternative="two-sided")
    with open(args.summary, "w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"n_manifest\t{len(samples)}\n")
        handle.write(f"n_result_rows\t{len(rows)}\n")
        handle.write(f"n_called\t{len(called)}\n")
        handle.write(f"n_missing\t{len(rows) - len(called)}\n")
        handle.write(f"case_DEL\t{case_del}\ncase_REF\t{case_ref}\n")
        handle.write(f"control_DEL\t{control_del}\ncontrol_REF\t{control_ref}\n")
        handle.write(f"fisher_odds_ratio\t{odds_ratio}\n")
        handle.write(f"fisher_two_sided_p\t{p_value}\n")


if __name__ == "__main__":
    main()
