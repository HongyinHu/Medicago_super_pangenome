#!/usr/bin/env bash
set -uo pipefail

PAN="path/to/project/N_1.EDTA_single/03_panEDTA"
cd "$PAN" || exit 2
mkdir -p logs/structural_resume/reconciled
LOG="logs/run_panEDTA.nohup.log"
RLOG="logs/run_panEDTA.resume_after_homology.log"
MLOG="logs/run_panEDTA.structural_monitor.log"
echo $$ > panEDTA.pid
total="$(wc -l < genome.list)"
echo "STRUCTURAL_RESUME_MONITOR_RESTART $(date '+%F %T %Z') total=$total pid=$$" | tee -a "$LOG" "$RLOG" "$MLOG" >/dev/null

while true; do
  bash run_panEDTA_structural_reconcile.sh >/dev/null 2>&1 || true
  done_count="$(find logs/structural_resume -name '*.done' 2>/dev/null | sed 's#.*/##;s/\.done$//' | sort -u | wc -l)"
  done_names="$(find logs/structural_resume -name '*.done' 2>/dev/null | sed 's#.*/##;s/\.done$//' | sort -u)"
  fail_count=0
  for ff in $(find logs/structural_resume -name '*.failed' 2>/dev/null); do
    bn="$(basename "$ff" .failed)"
    printf '%s\n' "$done_names" | grep -qx "$bn" && continue
    fail_count=$((fail_count+1))
  done
  active_count="$(ps -u $USER -ww -o cmd | grep -E 'EDTA.pl .*--step final' | grep -v grep | wc -l)"
  msg="$(date '+%F %T') structural_done=$done_count/$total failed_markers=$fail_count active_final=$active_count"
  echo "$msg" | tee -a "$LOG" "$RLOG" "$MLOG" >/dev/null
  if [ "$done_count" -ge "$total" ]; then
    echo "STRUCTURAL_RESUME_COMPLETE $(date '+%F %T %Z')" | tee -a "$LOG" "$RLOG" "$MLOG" >/dev/null
    break
  fi
  sleep 600
done