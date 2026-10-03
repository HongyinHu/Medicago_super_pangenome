#!/usr/bin/env bash
set -euo pipefail

OUT=path/to/project/N_3.call_SV/03.integrated_panSV
mkdir -p "$OUT/logs" "$OUT/results" "$OUT/summary"
LOG="$OUT/logs/run_integrate_panSV.$(date +%Y%m%d_%H%M%S).log"

{
  echo "HOST=$(hostname) DATE=$(date '+%F %T %Z')"
  echo "OUT=$OUT"
  python3 -u "$OUT/scripts/integrate_panSV.py"
} 2>&1 | tee "$LOG"
