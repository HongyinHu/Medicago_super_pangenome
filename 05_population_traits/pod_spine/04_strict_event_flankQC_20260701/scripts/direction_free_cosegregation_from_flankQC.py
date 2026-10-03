#!/usr/bin/env python3
"""Direction-free cosegregation summary from per-sample flank-QC SV calls."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


DEFAULT_RUN = Path(
    "path/to/project/"
    "N_4.pod_spiny/04_strict_event_flankQC_20260701"
)

MODELS = {
    "spiny_present": {"1": "present", "0": "absent"},
    "spiny_absent": {"1": "absent", "0": "present"},
}
VALID_CALLS = {"present", "absent"}


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def reference_self_sample(ref: str) -> str:
    if ref == "Msa":
        return "genome_Msa"
    if ref == "R108":
        return "genome_R108"
    return f"genome_{ref}"


def evaluate_model(rows: list[dict[str, str]], model: str) -> dict[str, object]:
    expected_by_pheno = MODELS[model]
    ok = 0
    conflict = 0
    ambiguous = 0
    na = 0
    total = 0
    spiny_ok = spiny_total = 0
    spineless_ok = spineless_total = 0
    conflict_samples = []
    ambiguous_samples = []
    for row in rows:
        pheno = row["spiny_primary"]
        if pheno not in expected_by_pheno:
            continue
        total += 1
        expected = expected_by_pheno[pheno]
        call = row["strict_event_call"]
        sample = row["sample"]
        if pheno == "1":
            spiny_total += 1
        else:
            spineless_total += 1
        if call == expected:
            ok += 1
            if pheno == "1":
                spiny_ok += 1
            else:
                spineless_ok += 1
        elif call in VALID_CALLS:
            conflict += 1
            conflict_samples.append(f"{sample}:{call}->expected_{expected}")
        else:
            if call == "NA":
                na += 1
            else:
                ambiguous += 1
            ambiguous_samples.append(f"{sample}:{call}->expected_{expected}")
    return {
        "model": model,
        "ok": ok,
        "conflict": conflict,
        "ambiguous": ambiguous,
        "na": na,
        "total": total,
        "spiny_ok": spiny_ok,
        "spiny_total": spiny_total,
        "spineless_ok": spineless_ok,
        "spineless_total": spineless_total,
        "conflict_samples": ",".join(conflict_samples) if conflict_samples else "NA",
        "ambiguous_samples": ",".join(ambiguous_samples) if ambiguous_samples else "NA",
    }


def choose_model(rows: list[dict[str, str]]) -> dict[str, object]:
    scores = [evaluate_model(rows, model) for model in MODELS]
    scores.sort(
        key=lambda item: (
            int(item["conflict"]),
            int(item["ambiguous"]) + int(item["na"]),
            -int(item["ok"]),
            item["model"],
        )
    )
    return scores[0]


def verdict(score: dict[str, object]) -> str:
    if int(score["total"]) == 0:
        return "no_core_samples"
    if int(score["spiny_total"]) == 0 or int(score["spineless_total"]) == 0:
        return "missing_one_phenotype_group"
    if int(score["conflict"]) == 0 and int(score["ambiguous"]) == 0 and int(score["na"]) == 0:
        return "direction_free_strict_pass"
    if int(score["conflict"]) == 0:
        return "direction_free_manual_missing_or_ambiguous"
    return "direction_free_fail_conflict"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    parser.add_argument(
        "--exclude-ref-self",
        action="store_true",
        default=True,
        help="Exclude the sample that is the current reference from core scoring.",
    )
    args = parser.parse_args()

    run = args.run_dir
    rows = load_tsv(run / "results" / "per_sample_event_flankQC.tsv")
    by_candidate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_candidate[row["candidate_id"]].append(row)

    out_rows = []
    per_model_rows = []
    for cid, cand_rows in sorted(by_candidate.items()):
        first = cand_rows[0]
        ref_self = reference_self_sample(first["ref"])
        core_rows = []
        sensitivity_calls = []
        for row in cand_rows:
            sample_call = f"{row['sample']}:{row['strict_event_call']}"
            if row.get("use_core_screen") == "1" and not (
                args.exclude_ref_self and row["sample"] == ref_self
            ):
                core_rows.append(row)
            else:
                sensitivity_calls.append(sample_call)
        model_scores = [evaluate_model(core_rows, model) for model in MODELS]
        for score in model_scores:
            per_model_rows.append(
                {
                    "candidate_id": cid,
                    "MetaSV_ID": first["MetaSV_ID"],
                    "ref": first["ref"],
                    "locus": f"{first['chrom']}:{first['start']}-{first['end']}",
                    "SVTYPE": first["SVTYPE"],
                    "SVLEN": first["SVLEN"],
                    "best_gene_id": first.get("best_gene_id", ""),
                    **{k: str(v) for k, v in score.items()},
                }
            )
        best = choose_model(core_rows)
        out_rows.append(
            {
                "candidate_id": cid,
                "MetaSV_ID": first["MetaSV_ID"],
                "ref": first["ref"],
                "locus": f"{first['chrom']}:{first['start']}-{first['end']}",
                "SVTYPE": first["SVTYPE"],
                "SVLEN": first["SVLEN"],
                "best_gene_id": first.get("best_gene_id", ""),
                "original_association_direction": first.get("association_direction", ""),
                "inferred_model": str(best["model"]),
                "inferred_spiny_state": MODELS[str(best["model"])]["1"],
                "inferred_spineless_state": MODELS[str(best["model"])]["0"],
                "core_ok": str(best["ok"]),
                "core_total": str(best["total"]),
                "core_conflict": str(best["conflict"]),
                "core_ambiguous": str(best["ambiguous"]),
                "core_NA": str(best["na"]),
                "spiny_ok": f"{best['spiny_ok']}/{best['spiny_total']}",
                "spineless_ok": f"{best['spineless_ok']}/{best['spineless_total']}",
                "conflict_samples": str(best["conflict_samples"]),
                "ambiguous_samples": str(best["ambiguous_samples"]),
                "sensitivity_calls": ";".join(sensitivity_calls) if sensitivity_calls else "NA",
                "direction_free_verdict": verdict(best),
            }
        )

    summary = run / "results" / "candidate_direction_free_summary.tsv"
    per_model = run / "results" / "candidate_direction_free_per_model.tsv"
    strict = run / "results" / "candidate_direction_free_strict_pass.tsv"
    counts = run / "summary" / "direction_free_counts.tsv"

    fieldnames = list(out_rows[0].keys())
    with summary.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    model_fieldnames = list(per_model_rows[0].keys())
    with per_model.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=model_fieldnames)
        writer.writeheader()
        writer.writerows(per_model_rows)

    with strict.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows([r for r in out_rows if r["direction_free_verdict"] == "direction_free_strict_pass"])

    count_by_verdict = defaultdict(int)
    count_by_model = defaultdict(int)
    for row in out_rows:
        count_by_verdict[row["direction_free_verdict"]] += 1
        count_by_model[row["inferred_model"]] += 1
    with counts.open("w") as out:
        out.write("category\tname\tcount\n")
        for name, count in sorted(count_by_verdict.items()):
            out.write(f"verdict\t{name}\t{count}\n")
        for name, count in sorted(count_by_model.items()):
            out.write(f"model\t{name}\t{count}\n")

    print(f"summary={summary}")
    print(f"strict_pass={strict}")
    print(f"counts={counts}")


if __name__ == "__main__":
    main()
