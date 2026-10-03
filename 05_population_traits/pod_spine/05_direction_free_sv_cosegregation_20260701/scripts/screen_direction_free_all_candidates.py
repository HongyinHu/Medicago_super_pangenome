#!/usr/bin/env python3
"""Direction-free SV/phenotype cosegregation screen for pod-spiny candidates."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from collections import Counter
from pathlib import Path


DEFAULT_N4 = Path("path/to/project/N_4.pod_spiny")
DEFAULT_OUT = DEFAULT_N4 / "05_direction_free_sv_cosegregation_20260701"
DEFAULT_INPUT = (
    DEFAULT_N4
    / "results"
    / "meta_publication"
    / "validation_priority_candidates.with_coords.genes.te.tsv"
)
DEFAULT_SAMPLE_CONFIG = (
    DEFAULT_N4 / "04_strict_event_flankQC_20260701" / "config" / "samples_extended.tsv"
)

MODELS = {
    "spiny_present": {"1": "present", "0": "absent"},
    "spiny_absent": {"1": "absent", "0": "present"},
}
VALID_CALLS = {"present", "absent"}


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def norm(value: str | None) -> str:
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value in {".", "NA", "na", "None", "none"} else value


def split_samples(value: str | None) -> set[str]:
    value = norm(value)
    if not value:
        return set()
    return {x.strip() for x in value.split(",") if x.strip()}


def as_int(value: str | None, default: int = 0) -> int:
    value = norm(value)
    if not value:
        return default
    return int(float(value))


def as_float(value: str | None, default: float = 0.0) -> float:
    value = norm(value)
    if not value:
        return default
    return float(value)


def choose_coord(row: dict[str, str]) -> tuple[str, str, int, int, str, str]:
    """Return igv_ref, chrom, start, end, candidate_locus, refs_with_coords."""
    refs = []
    if norm(row.get("RefA_CHROM")):
        refs.append("Msa")
    if norm(row.get("RefB_CHROM")):
        refs.append("R108")
    best_ref = norm(row.get("best_gene_ref"))
    if best_ref == "RefA" and "Msa" in refs:
        ref = "Msa"
    elif best_ref == "RefB" and "R108" in refs:
        ref = "R108"
    elif "Msa" in refs:
        ref = "Msa"
    elif "R108" in refs:
        ref = "R108"
    else:
        return "", "", 0, 0, "", ""
    if ref == "Msa":
        chrom = row["RefA_CHROM"]
        start = as_int(row["RefA_POS"])
        end = as_int(row["RefA_END"])
        locus = f"RefA:{chrom}:{start}-{end}"
    else:
        chrom = row["RefB_CHROM"]
        start = as_int(row["RefB_POS"])
        end = as_int(row["RefB_END"])
        locus = f"RefB:{chrom}:{start}-{end}"
    return ref, chrom, start, end, locus, ",".join(refs)


def sample_calls_from_row(row: dict[str, str]) -> dict[str, str]:
    calls: dict[str, str] = {}
    for sample in split_samples(row.get("spiny_present_samples")) | split_samples(row.get("spineless_present_samples")):
        calls[sample] = "present"
    for sample in split_samples(row.get("spiny_absent_samples")) | split_samples(row.get("spineless_absent_samples")):
        calls[sample] = "absent"
    return calls


def references_to_exclude(refs_with_coords: str) -> set[str]:
    refs = set(refs_with_coords.split(",")) if refs_with_coords else set()
    excluded = set()
    if "Msa" in refs:
        excluded.add("genome_Msa")
    if "R108" in refs:
        excluded.add("genome_R108")
    return excluded


def evaluate_model(
    row: dict[str, str],
    sample_rows: list[dict[str, str]],
    model: str,
    refs_with_coords: str,
) -> dict[str, object]:
    calls = sample_calls_from_row(row)
    excluded_refs = references_to_exclude(refs_with_coords)
    expected_by_pheno = MODELS[model]
    ok = conflict = missing = ambiguous = 0
    spiny_ok = spiny_total = 0
    spineless_ok = spineless_total = 0
    conflict_samples = []
    missing_samples = []
    used_samples = []
    excluded_samples = []
    for sample_row in sample_rows:
        sample = sample_row["sample"]
        pheno = sample_row["spiny_primary"]
        if sample_row.get("use_core_screen") != "1":
            continue
        if sample in excluded_refs:
            excluded_samples.append(sample)
            continue
        if pheno not in expected_by_pheno:
            continue
        expected = expected_by_pheno[pheno]
        call = calls.get(sample, "NA")
        used_samples.append(sample)
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
        elif call == "NA":
            missing += 1
            missing_samples.append(f"{sample}:NA->expected_{expected}")
        else:
            ambiguous += 1
            missing_samples.append(f"{sample}:{call}->expected_{expected}")
    return {
        "model": model,
        "ok": ok,
        "conflict": conflict,
        "missing": missing,
        "ambiguous": ambiguous,
        "total": len(used_samples),
        "spiny_ok": spiny_ok,
        "spiny_total": spiny_total,
        "spineless_ok": spineless_ok,
        "spineless_total": spineless_total,
        "conflict_samples": ",".join(conflict_samples) if conflict_samples else "NA",
        "missing_samples": ",".join(missing_samples) if missing_samples else "NA",
        "excluded_reference_samples": ",".join(sorted(set(excluded_samples))) if excluded_samples else "NA",
    }


def choose_model(scores: list[dict[str, object]]) -> dict[str, object]:
    return sorted(
        scores,
        key=lambda s: (
            int(s["conflict"]),
            int(s["missing"]) + int(s["ambiguous"]),
            -int(s["ok"]),
            s["model"],
        ),
    )[0]


def verdict(score: dict[str, object]) -> str:
    if int(score["spiny_total"]) < 3 or int(score["spineless_total"]) < 8:
        return "insufficient_core_samples"
    if int(score["conflict"]) == 0 and int(score["missing"]) == 0 and int(score["ambiguous"]) == 0:
        return "direction_free_exact_core"
    if int(score["conflict"]) == 0 and int(score["missing"]) + int(score["ambiguous"]) <= 2:
        return "direction_free_near_core_missing_le2"
    if int(score["conflict"]) == 0:
        return "direction_free_manual_missing"
    return "direction_free_fail_conflict"


def igv_priority(row: dict[str, object]) -> float:
    score = as_float(str(row.get("final_validation_score", "")), 0.0)
    if row.get("direction_free_verdict") == "direction_free_exact_core":
        score += 1000
    elif row.get("direction_free_verdict") == "direction_free_near_core_missing_le2":
        score += 500
    if str(row.get("SVTYPE")) == "DEL":
        score += 50
    if str(row.get("best_gene_context")) == "overlap_gene":
        score += 30
    if str(row.get("TE_overlap_flag")) == "no_TE_overlap":
        score += 20
    if 50 <= abs(as_int(str(row.get("SVLEN")), 0)) <= 5000:
        score += 10
    return score


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--sample-config", type=Path, default=DEFAULT_SAMPLE_CONFIG)
    parser.add_argument("--max-igv", type=int, default=80)
    parser.add_argument("--igv-verdicts", default="direction_free_exact_core,direction_free_near_core_missing_le2")
    args = parser.parse_args()

    out_dir = args.out_dir
    for sub in ("config", "results", "summary", "scripts", "logs"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.sample_config, out_dir / "config" / "samples_extended.tsv")
    shutil.copy2(args.input, out_dir / "config" / "input_candidates.with_coords.genes.te.tsv")

    sample_rows = load_tsv(args.sample_config)
    rows = load_tsv(args.input)
    out_rows = []
    per_model_rows = []
    for idx, row in enumerate(rows, 1):
        igv_ref, chrom, start, end, locus, refs_with_coords = choose_coord(row)
        if not igv_ref:
            continue
        cid = f"{int(row.get('final_validation_rank') or idx):03d}_{row['MetaSV_ID']}"
        scores = [evaluate_model(row, sample_rows, model, refs_with_coords) for model in MODELS]
        for score in scores:
            per_model_rows.append(
                {
                    "candidate_id": cid,
                    "MetaSV_ID": row["MetaSV_ID"],
                    "igv_ref": igv_ref,
                    "locus": locus,
                    "SVTYPE": row["SVTYPE"],
                    "SVLEN": row["SVLEN"],
                    **{k: str(v) for k, v in score.items()},
                }
            )
        best = choose_model(scores)
        merged = {
            **row,
            "candidate_id": cid,
            "igv_ref": igv_ref,
            "igv_chrom": chrom,
            "igv_start": start,
            "igv_end": end,
            "igv_locus": locus,
            "refs_with_coords": refs_with_coords,
            "inferred_model": best["model"],
            "inferred_spiny_state": MODELS[str(best["model"])]["1"],
            "inferred_spineless_state": MODELS[str(best["model"])]["0"],
            "core_ok": best["ok"],
            "core_total": best["total"],
            "core_conflict": best["conflict"],
            "core_missing": best["missing"],
            "core_ambiguous": best["ambiguous"],
            "spiny_ok": f"{best['spiny_ok']}/{best['spiny_total']}",
            "spineless_ok": f"{best['spineless_ok']}/{best['spineless_total']}",
            "conflict_samples": best["conflict_samples"],
            "missing_samples": best["missing_samples"],
            "excluded_reference_samples": best["excluded_reference_samples"],
            "direction_free_verdict": verdict(best),
        }
        merged["igv_priority_score"] = f"{igv_priority(merged):.3f}"
        out_rows.append(merged)

    keep_fields = [
        "candidate_id",
        "MetaSV_ID",
        "igv_ref",
        "igv_chrom",
        "igv_start",
        "igv_end",
        "igv_locus",
        "refs_with_coords",
        "SVTYPE",
        "SVLEN",
        "best_gene_ref",
        "best_gene_context",
        "best_gene_id",
        "best_gene_locus",
        "TE_overlap_flag",
        "best_TE_overlap_fraction",
        "best_TE_classes",
        "final_validation_rank",
        "final_validation_score",
        "priority_class",
        "validation_bucket",
        "tier",
        "support_models",
        "n_support_models",
        "original_association_direction",
        "association_direction",
        "inferred_model",
        "inferred_spiny_state",
        "inferred_spineless_state",
        "core_ok",
        "core_total",
        "core_conflict",
        "core_missing",
        "core_ambiguous",
        "spiny_ok",
        "spineless_ok",
        "conflict_samples",
        "missing_samples",
        "excluded_reference_samples",
        "direction_free_verdict",
        "igv_priority_score",
        "candidate_locus",
        "RefA_CHROM",
        "RefA_POS",
        "RefA_END",
        "RefB_CHROM",
        "RefB_POS",
        "RefB_END",
        "RefA_SV_ID",
        "RefB_SV_ID",
        "match_status",
        "match_method",
        "confidence",
        "spiny_present_samples",
        "spineless_present_samples",
        "spiny_absent_samples",
        "spineless_absent_samples",
    ]
    # Preserve only fields that exist or are generated.
    keep_fields = [f for f in keep_fields if any(f in row for row in out_rows)]
    out_rows.sort(
        key=lambda r: (
            {
                "direction_free_exact_core": 0,
                "direction_free_near_core_missing_le2": 1,
                "direction_free_manual_missing": 2,
                "direction_free_fail_conflict": 3,
            }.get(str(r["direction_free_verdict"]), 9),
            -as_float(str(r.get("igv_priority_score")), 0.0),
            as_int(str(r.get("final_validation_rank")), 10**9),
        )
    )

    summary_path = out_dir / "results" / "direction_free_all_candidates.tsv"
    write_tsv(summary_path, out_rows, keep_fields)
    exact_rows = [r for r in out_rows if r["direction_free_verdict"] == "direction_free_exact_core"]
    near_rows = [
        r
        for r in out_rows
        if r["direction_free_verdict"]
        in {"direction_free_exact_core", "direction_free_near_core_missing_le2"}
    ]
    write_tsv(out_dir / "results" / "direction_free_exact_core.tsv", exact_rows, keep_fields)
    write_tsv(out_dir / "results" / "direction_free_near_or_exact.tsv", near_rows, keep_fields)
    write_tsv(
        out_dir / "results" / "direction_free_per_model.tsv",
        per_model_rows,
        list(per_model_rows[0].keys()) if per_model_rows else [],
    )

    igv_verdicts = {x.strip() for x in args.igv_verdicts.split(",") if x.strip()}
    igv_candidates = [
        r
        for r in out_rows
        if r["direction_free_verdict"] in igv_verdicts
        and r.get("igv_ref")
        and r.get("igv_chrom")
        and r.get("SVTYPE") == "DEL"
        and 50 <= abs(as_int(str(r.get("SVLEN")), 0)) <= 5000
    ]
    igv_candidates.sort(
        key=lambda r: (
            -as_float(str(r.get("igv_priority_score")), 0.0),
            as_int(str(r.get("final_validation_rank")), 10**9),
        )
    )
    if args.max_igv > 0:
        igv_candidates = igv_candidates[: args.max_igv]
    write_tsv(out_dir / "results" / "direction_free_igv_selected.tsv", igv_candidates, keep_fields)

    # Compatibility tables for the IGV report builder used in step 04.
    config_rows = []
    summary_rows = []
    for r in igv_candidates:
        ref = "Msa" if r["igv_ref"] == "Msa" else "R108"
        chrom, start, end = r["igv_chrom"], int(r["igv_start"]), int(r["igv_end"])
        config_rows.append(
            {
                "candidate_id": r["candidate_id"],
                "MetaSV_ID": r["MetaSV_ID"],
                "ref": ref,
                "chrom": chrom,
                "start": start,
                "end": end,
                "view_region": f"{chrom}:{max(1, start - 10000)}-{end + 10000}",
                "SVTYPE": r["SVTYPE"],
                "SVLEN": r["SVLEN"],
                "association_direction": r["inferred_model"],
                "best_gene_id": r.get("best_gene_id", ""),
                "candidate_locus": r["igv_locus"],
            }
        )
        summary_rows.append(
            {
                "candidate_id": r["candidate_id"],
                "MetaSV_ID": r["MetaSV_ID"],
                "ref": ref,
                "locus": f"{chrom}:{start}-{end}",
                "SVTYPE": r["SVTYPE"],
                "SVLEN": r["SVLEN"],
                "association_direction": r["inferred_model"],
                "best_gene_id": r.get("best_gene_id", ""),
                "core_expected_present_ok": r.get("spiny_ok", ""),
                "core_expected_absent_ok": r.get("spineless_ok", ""),
                "core_present_calls": "",
                "core_absent_calls": "",
                "core_ambiguous_calls": r.get("core_ambiguous", ""),
                "core_NA_calls": r.get("core_missing", ""),
                "core_conflict_samples": r.get("conflict_samples", ""),
                "core_unknown_samples": r.get("missing_samples", ""),
                "core_low_flank_samples": "NA",
                "sensitivity_calls": f"excluded_reference_samples={r.get('excluded_reference_samples', 'NA')}",
                "strict_flankQC_verdict": r["direction_free_verdict"],
            }
        )
    write_tsv(
        out_dir / "config" / "candidates_29.tsv",
        config_rows,
        [
            "candidate_id",
            "MetaSV_ID",
            "ref",
            "chrom",
            "start",
            "end",
            "view_region",
            "SVTYPE",
            "SVLEN",
            "association_direction",
            "best_gene_id",
            "candidate_locus",
        ],
    )
    write_tsv(
        out_dir / "results" / "candidate_strict_flankQC_summary.tsv",
        summary_rows,
        [
            "candidate_id",
            "MetaSV_ID",
            "ref",
            "locus",
            "SVTYPE",
            "SVLEN",
            "association_direction",
            "best_gene_id",
            "core_expected_present_ok",
            "core_expected_absent_ok",
            "core_present_calls",
            "core_absent_calls",
            "core_ambiguous_calls",
            "core_NA_calls",
            "core_conflict_samples",
            "core_unknown_samples",
            "core_low_flank_samples",
            "sensitivity_calls",
            "strict_flankQC_verdict",
        ],
    )

    counts = Counter(r["direction_free_verdict"] for r in out_rows)
    model_counts = Counter(r["inferred_model"] for r in out_rows)
    type_counts = Counter(r["SVTYPE"] for r in out_rows)
    with (out_dir / "summary" / "direction_free_counts.tsv").open("w") as out:
        out.write("category\tname\tcount\n")
        for name, count in sorted(counts.items()):
            out.write(f"verdict\t{name}\t{count}\n")
        for name, count in sorted(model_counts.items()):
            out.write(f"model\t{name}\t{count}\n")
        for name, count in sorted(type_counts.items()):
            out.write(f"SVTYPE\t{name}\t{count}\n")
        out.write(f"igv_selected\tDEL_50bp_5kb_near_or_exact\t{len(igv_candidates)}\n")
    with (out_dir / "README.md").open("w") as out:
        out.write("# 05 direction-free SV phenotype cosegregation\n\n")
        out.write("Input: meta_publication validation-priority candidates with coordinates and annotation.\n")
        out.write("Method: test both spiny=present and spiny=absent models, choose the model with fewer conflicts, then fewer missing calls.\n")
        out.write("Reference-self samples are excluded from core scoring when the event has that reference coordinate.\n")
        out.write("IGV selection is limited to DEL 50bp-5kb among exact/near direction-free candidates.\n")
    print(f"input_rows={len(rows)}")
    print(f"scored_rows={len(out_rows)}")
    print(f"exact_core={counts.get('direction_free_exact_core', 0)}")
    print(f"near_or_exact={len(near_rows)}")
    print(f"igv_selected={len(igv_candidates)}")
    print(f"summary={summary_path}")
    print(f"igv_table={out_dir / 'results' / 'direction_free_igv_selected.tsv'}")


if __name__ == "__main__":
    main()
