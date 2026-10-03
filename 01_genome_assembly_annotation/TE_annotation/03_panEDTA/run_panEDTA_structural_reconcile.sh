#!/usr/bin/env bash
set -uo pipefail
PAN="path/to/project/N_1.EDTA_single/03_panEDTA"
cd "$PAN" || exit 2
mkdir -p logs/structural_resume/reconciled
while read -r g rest; do
  b="$(basename "$g")"
  [ -z "$b" ] && continue
  sum="$b.mod.EDTA.TEanno.sum"
  latest_log="$(ls -t logs/structural_resume/*/$b.EDTA_final.20260623_*.log 2>/dev/null | head -1)"
  if [ -n "$latest_log" ] && [ -s "$sum" ] && grep -q 'Evaluation of TE annotation finished' "$latest_log"; then
    rm -f logs/structural_resume/*/$b.failed 2>/dev/null || true
    echo "done $(date '+%F %T') sum=$sum log=$latest_log" > "logs/structural_resume/reconciled/$b.done"
  fi
done < genome.list