#!/usr/bin/env bash
set -euo pipefail

OUT=path/to/project/N_3.call_SV/04.integrated_panSV_Msa_single_ref
mkdir -p "$OUT/logs" "$OUT/results" "$OUT/summary"
LOG="$OUT/logs/run_integrate_panSV_Msa_single_ref.$(date +%Y%m%d_%H%M%S).log"

{
  echo "HOST=$(hostname) DATE=$(date '+%F %T %Z')"
  echo "OUT=$OUT"
  python3 -u "$OUT/scripts/integrate_panSV_Msa_single_ref.py"
} 2>&1 | tee "$LOG"
