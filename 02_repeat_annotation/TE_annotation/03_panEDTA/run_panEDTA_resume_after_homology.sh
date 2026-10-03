#!/usr/bin/env bash
set -uo pipefail

PAN="path/to/project/N_1.EDTA_single/03_panEDTA"
LOG="$PAN/logs/run_panEDTA.resume_after_homology.log"
cd "$PAN" || exit 2
echo $$ > panEDTA_resume_waiter.pid

total="$(wc -l < genome.list)"
echo "START_WAIT $(date '+%F %T %Z') host=$(hostname) total=$total" >> "$LOG"

while true; do
  out_count="$(find . -maxdepth 1 -name '*.mod.panEDTA.out' | wc -l)"
  lock_count="$(find logs/homology_locks -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)"
  echo "$(date '+%F %T') homology_out=$out_count/$total locks=$lock_count" >> "$LOG"
  if [ "$out_count" -ge "$total" ] && [ "$lock_count" -eq 0 ]; then
    break
  fi
  sleep 600
done

if [ -s panEDTA.pid ]; then
  oldpid="$(cat panEDTA.pid)"
  if ps -p "$oldpid" >/dev/null 2>&1; then
    echo "$(date '+%F %T') panEDTA already running pid=$oldpid; waiter exits" >> "$LOG"
    exit 0
  fi
fi

echo "$(date '+%F %T') launching run_panEDTA.sh after homology completion" >> "$LOG"
echo $$ > panEDTA.pid
export THREADS="${THREADS:-20}"
exec bash run_panEDTA.sh >> "$LOG" 2>&1
