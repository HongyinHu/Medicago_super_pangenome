#!/usr/bin/env bash
set -euo pipefail

OUT="path/to/project/N_4.pod_spiny/05_direction_free_sv_cosegregation_20260701"
N4="path/to/project/N_4.pod_spiny"
MAIN="path/to/project/N_3.call_SV/01.read_based_dualref_hifi"
PY="path/to/home/anaconda3/bin/python3"

mkdir -p "$OUT/logs"

"$PY" "$OUT/scripts/make_igv_reports_direction_free.py" \
  --run-dir "$OUT" \
  --n4-dir "$N4" \
  --main-sv-dir "$MAIN" \
  --only-ref ALL \
  --verdicts strict_pass,needs_manual_review_missing_or_ambiguous \
  --jobs 16 \
  --report-jobs 4 \
  > "$OUT/logs/make_igv_reports_direction_free.out" \
  2> "$OUT/logs/make_igv_reports_direction_free.err"
