#!/usr/bin/env python3
from pathlib import Path
import csv

run = Path("path/to/project/N_4.pod_spiny/05_direction_free_sv_cosegregation_20260701")
files = [
    ("all", run / "results/direction_free_all_candidates.tsv"),
    ("igv_selected", run / "results/direction_free_igv_selected.tsv"),
    ("flank_summary", run / "results/candidate_strict_flankQC_summary.tsv"),
    ("strict", run / "results/candidate_strict_pass.tsv"),
    ("manual", run / "results/candidate_manual_review.tsv"),
    ("fail", run / "results/candidate_fail_conflict.tsv"),
    ("manifest", run / "igv_reports_html/igv_reports_manifest.tsv"),
]
keep = [
    "candidate_id",
    "MetaSV_ID",
    "igv_ref",
    "ref",
    "igv_locus",
    "locus",
    "SVTYPE",
    "SVLEN",
    "best_gene_id",
    "best_gene_context",
    "strict_flankQC_verdict",
    "association_direction",
    "direction_free_verdict",
    "core_conflict",
    "core_missing",
    "core_ok",
    "core_total",
]

for label, path in files:
    print(f"## {label}\t{path}")
    if not path.exists():
        print("MISSING")
        continue
    hits = 0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            text = "\t".join(str(v) for v in row.values())
            if "Chr23997" in text or "23997" in text:
                hits += 1
                print("\t".join(f"{k}={row.get(k, '')}" for k in keep if k in row))
    print(f"hits\t{hits}")
