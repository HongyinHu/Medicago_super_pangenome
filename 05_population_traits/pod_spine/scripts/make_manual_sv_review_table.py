#!/usr/bin/env python3
"""Create strict manual-review tables for pod-spiny IGV SV candidates."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict, Counter
from pathlib import Path

SAMPLES = [
    "genome_410", "genome_474", "genome_Mpo", "genome_R108",
    "genome_395", "genome_454", "genome_461", "genome_468", "genome_472",
    "genome_482", "genome_M22", "genome_Mar", "genome_Mru", "genome_Msa", "genome_ZM4",
]
SPINY = {"genome_410", "genome_474", "genome_Mpo", "genome_R108"}
SPINELESS = set(SAMPLES) - SPINY


def fnum(row, key, default=0.0):
    try:
        return float(row.get(key, default) or default)
    except Exception:
        return default


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def classify_sample(row: dict[str, str]) -> tuple[str, str]:
    """Strict exact-interval DEL support classification from read-support metrics."""
    reads = fnum(row, "reads_in_window")
    del_reads = fnum(row, "del_cigar_support_reads")
    softclip = fnum(row, "softclip_near_breakpoint_reads")
    ratio = fnum(row, "event_to_flank_cov_ratio", 999.0)
    flank = fnum(row, "flank_cov_mean")

    if reads < 8 or flank < 5:
        return "ambiguous", "low_read_depth"

    # Clear deletion evidence in the exact candidate interval.
    if del_reads >= 3 and ratio <= 0.55:
        return "support", f"DEL_reads={del_reads:g};cov_ratio={ratio:.3g}"
    if del_reads >= 8 and ratio <= 0.75:
        return "support", f"many_DEL_reads={del_reads:g};cov_ratio={ratio:.3g}"
    if ratio <= 0.25 and (del_reads >= 1 or softclip >= 2):
        return "support", f"coverage_drop;DEL_reads={del_reads:g};softclip={softclip:g};cov_ratio={ratio:.3g}"

    # Strong reference-like evidence. Keep the threshold conservative.
    if del_reads <= 1 and ratio >= 0.75:
        return "no_support", f"reference_like;DEL_reads={del_reads:g};cov_ratio={ratio:.3g}"

    return "ambiguous", f"borderline;DEL_reads={del_reads:g};softclip={softclip:g};cov_ratio={ratio:.3g}"


def candidate_decision(spiny_calls: Counter, spineless_calls: Counter, context: str) -> tuple[str, str]:
    spiny_support = spiny_calls["support"]
    spiny_amb = spiny_calls["ambiguous"]
    spineless_support = spineless_calls["support"]
    spineless_amb = spineless_calls["ambiguous"]
    spineless_no = spineless_calls["no_support"]

    reasons = []
    if spineless_support > 0:
        reasons.append(f"spineless_exact_support={spineless_support}")
    if spineless_amb >= 3:
        reasons.append(f"spineless_many_ambiguous={spineless_amb}")
    if spiny_support < 3:
        reasons.append(f"spiny_support_low={spiny_support}/4")
    if context in {"CDS", "CDS+intron_boundary", "exon_nonCDS", "exon+intron_boundary"}:
        reasons.append(f"functional_context={context}")
    elif context == "intron":
        reasons.append("intronic")
    elif context in {"upstream", "downstream"}:
        reasons.append(context)

    if spiny_support >= 4 and spineless_support == 0 and spineless_amb <= 1:
        return "Tier1_exact_interval", ";".join(reasons or ["4/4_spiny;0/11_spineless"])
    if spiny_support >= 3 and spineless_support == 0 and spineless_amb <= 2:
        return "Tier2_manual_review", ";".join(reasons or ["3-4_spiny;0_spineless_support"])
    if spiny_support >= 3 and spineless_support <= 1 and spineless_no >= 8:
        return "Tier3_secondary_only", ";".join(reasons or ["partial_pattern"])
    return "Exclude_or_low_priority", ";".join(reasons or ["weak_or_nonsegregating_pattern"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--package-dir", type=Path, default=Path("path/to/project/N_4.pod_spiny/02_IGV_validation_29SV"))
    args = ap.parse_args()
    pkg = args.package_dir
    metrics = load_tsv(pkg / "read_support_metrics.tsv")
    candidate_rows = {r["candidate_id"]: r for r in load_tsv(pkg / "config" / "candidates_29.tsv")}
    gene_context_path = pkg / "igv_reports_html" / "candidate_gene_model_context.tsv"
    gene_context = {r["candidate_id"]: r for r in load_tsv(gene_context_path)} if gene_context_path.exists() else {}
    support_summary_path = pkg / "candidate_read_support_summary.tsv"
    support_summary = {r["candidate_id"]: r for r in load_tsv(support_summary_path)} if support_summary_path.exists() else {}

    rows_by_candidate = defaultdict(list)
    for row in metrics:
        rows_by_candidate[row["candidate_id"]].append(row)

    out_dir = pkg / "manual_review_strict"
    out_dir.mkdir(parents=True, exist_ok=True)

    per_sample_path = out_dir / "strict_per_sample_calls.tsv"
    summary_path = out_dir / "strict_candidate_decision_summary.tsv"
    matrix_path = out_dir / "manual_review_matrix.tsv"
    shortlist_path = out_dir / "tier1_tier2_review_shortlist.tsv"

    per_sample_fields = [
        "candidate_id", "MetaSV_ID", "sample", "phenotype", "auto_call", "auto_reason",
        "reads_in_window", "del_cigar_support_reads", "softclip_near_breakpoint_reads",
        "event_to_flank_cov_ratio", "manual_call", "manual_note",
    ]
    with per_sample_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=per_sample_fields)
        writer.writeheader()
        for cid in sorted(rows_by_candidate, key=lambda x: int(x.split("_")[0]) if x.split("_")[0].isdigit() else 9999):
            sample_rows = {r["sample"]: r for r in rows_by_candidate[cid]}
            for sample in SAMPLES:
                r = sample_rows.get(sample)
                if r is None:
                    writer.writerow({"candidate_id": cid, "sample": sample, "auto_call": "missing", "auto_reason": "no_metric_row"})
                    continue
                call, reason = classify_sample(r)
                writer.writerow({
                    "candidate_id": cid,
                    "MetaSV_ID": r["MetaSV_ID"],
                    "sample": sample,
                    "phenotype": r["phenotype"],
                    "auto_call": call,
                    "auto_reason": reason,
                    "reads_in_window": r["reads_in_window"],
                    "del_cigar_support_reads": r["del_cigar_support_reads"],
                    "softclip_near_breakpoint_reads": r["softclip_near_breakpoint_reads"],
                    "event_to_flank_cov_ratio": r["event_to_flank_cov_ratio"],
                    "manual_call": "",
                    "manual_note": "",
                })

    summary_rows = []
    for cid in sorted(rows_by_candidate, key=lambda x: int(x.split("_")[0]) if x.split("_")[0].isdigit() else 9999):
        calls = {}
        reasons = {}
        for r in rows_by_candidate[cid]:
            call, reason = classify_sample(r)
            calls[r["sample"]] = call
            reasons[r["sample"]] = reason
        spiny_counter = Counter(calls.get(s, "missing") for s in SPINY)
        spineless_counter = Counter(calls.get(s, "missing") for s in SPINELESS)
        cand = candidate_rows.get(cid, {})
        ctx = gene_context.get(cid, {})
        supp = support_summary.get(cid, {})
        context = ctx.get("gene_model_context", "")
        decision, decision_reason = candidate_decision(spiny_counter, spineless_counter, context)
        html_name = ""
        if cand:
            html_name = f"../igv_reports_html/{cid}_{cand.get('SVTYPE','SV')}_len{cand.get('SVLEN','NA')}_{cand.get('chrom','chr')}_{cand.get('start','start')}_{cand.get('end','end')}.html"
        summary_rows.append({
            "candidate_id": cid,
            "MetaSV_ID": cand.get("MetaSV_ID", ""),
            "ref": cand.get("ref", ""),
            "locus": f"{cand.get('chrom','')}:{cand.get('start','')}-{cand.get('end','')}",
            "SVTYPE": cand.get("SVTYPE", ""),
            "SVLEN": cand.get("SVLEN", ""),
            "association_direction": cand.get("association_direction", ""),
            "best_gene_id": cand.get("best_gene_id", ""),
            "gene_model_context": context,
            "context_gene_id": ctx.get("context_gene_id", ""),
            "context_gene_locus": ctx.get("context_gene_locus", ""),
            "context_distance_bp": ctx.get("context_distance_bp", ""),
            "original_support_verdict": supp.get("support_verdict", ""),
            "auto_spiny_support": str(spiny_counter["support"]),
            "auto_spiny_ambiguous": str(spiny_counter["ambiguous"]),
            "auto_spiny_no_support": str(spiny_counter["no_support"]),
            "auto_spineless_support": str(spineless_counter["support"]),
            "auto_spineless_ambiguous": str(spineless_counter["ambiguous"]),
            "auto_spineless_no_support": str(spineless_counter["no_support"]),
            "auto_decision": decision,
            "auto_decision_reason": decision_reason,
            "manual_event_consistency": "",
            "manual_final_decision": "",
            "manual_notes": "",
            "html_report": html_name,
            **{f"auto_{s}": calls.get(s, "missing") for s in SAMPLES},
        })

    summary_fields = [
        "candidate_id", "MetaSV_ID", "ref", "locus", "SVTYPE", "SVLEN", "association_direction",
        "best_gene_id", "gene_model_context", "context_gene_id", "context_gene_locus", "context_distance_bp",
        "original_support_verdict", "auto_spiny_support", "auto_spiny_ambiguous", "auto_spiny_no_support",
        "auto_spineless_support", "auto_spineless_ambiguous", "auto_spineless_no_support",
        "auto_decision", "auto_decision_reason", "manual_event_consistency", "manual_final_decision",
        "manual_notes", "html_report",
    ] + [f"auto_{s}" for s in SAMPLES]
    with summary_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    with shortlist_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=summary_fields)
        writer.writeheader()
        for row in summary_rows:
            if row["auto_decision"] in {"Tier1_exact_interval", "Tier2_manual_review"}:
                writer.writerow(row)

    matrix_fields = ["candidate_id", "MetaSV_ID", "locus", "SVLEN", "gene_model_context", "auto_decision"] + SAMPLES + ["manual_final_decision", "manual_notes"]
    with matrix_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=matrix_fields)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({
                "candidate_id": row["candidate_id"],
                "MetaSV_ID": row["MetaSV_ID"],
                "locus": row["locus"],
                "SVLEN": row["SVLEN"],
                "gene_model_context": row["gene_model_context"],
                "auto_decision": row["auto_decision"],
                **{s: row.get(f"auto_{s}", "") for s in SAMPLES},
                "manual_final_decision": "",
                "manual_notes": "",
            })

    c = Counter(r["auto_decision"] for r in summary_rows)
    print("Wrote", summary_path)
    print("Wrote", per_sample_path)
    print("Wrote", matrix_path)
    print("Wrote", shortlist_path)
    for k, v in c.most_common():
        print(k, v)


if __name__ == "__main__":
    main()
