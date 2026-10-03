#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="path/to/project/N_4.pod_spiny/04_strict_event_flankQC_20260701"
N4="path/to/project/N_4.pod_spiny"
N3="path/to/project/N_3.call_SV/01.read_based_dualref_hifi"

mkdir -p "${RUN_DIR}"/{config,logs,results,scripts,summary}

path/to/home/anaconda3/bin/python3 "${RUN_DIR}/scripts/flank_qc_event_screen.py" \
  --candidates "${RUN_DIR}/config/candidates_29.tsv" \
  --samples "${RUN_DIR}/config/samples_extended.tsv" \
  --bam-root "${N3}/03_per_sample" \
  --vcf-root "${N3}/03_per_sample" \
  --out-dir "${RUN_DIR}/results" \
  --samtools path/to/home/anaconda3/bin/samtools \
  --tabix path/to/home/anaconda3/bin/tabix \
  --flank 1000 \
  --pad 300 \
  --threads 8 \
  --min-flank-depth 5 \
  --min-flank-covered 0.70 \
  --max-breakpoint-dist 200 \
  --min-reciprocal-overlap 0.50 \
  > "${RUN_DIR}/logs/flank_qc_event_screen.out" \
  2> "${RUN_DIR}/logs/flank_qc_event_screen.err"

path/to/home/anaconda3/bin/python3 - <<'PY'
import csv
from collections import Counter
from pathlib import Path
run = Path("path/to/project/N_4.pod_spiny/04_strict_event_flankQC_20260701")
summary = run / "results" / "candidate_strict_flankQC_summary.tsv"
counts = Counter()
rows = []
with summary.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):
        counts[r["strict_flankQC_verdict"]] += 1
        rows.append(r)
with (run / "summary" / "flankQC_counts.tsv").open("w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["verdict", "count"])
    for k, v in sorted(counts.items()):
        w.writerow([k, v])
with (run / "summary" / "flankQC_top_lines.tsv").open("w", newline="") as out:
    fields = ["candidate_id", "MetaSV_ID", "locus", "SVLEN", "association_direction", "core_expected_present_ok", "core_expected_absent_ok", "core_conflict_samples", "core_unknown_samples", "strict_flankQC_verdict"]
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows[:50]:
        w.writerow(r)
PY

touch "${RUN_DIR}/summary/flank_qc_test.done"
