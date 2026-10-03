#!/usr/bin/env bash
set -uo pipefail

PAN="path/to/project/N_1.EDTA_single/03_panEDTA"
LIST_FILE="${1:?list file required}"
NODE_LABEL="${2:?node label required}"
LOGDIR="$PAN/logs/direct_homology/$NODE_LABEL"
LOCKDIR="$PAN/logs/homology_locks"
LOG="$LOGDIR/monitor.log"

mkdir -p "$LOGDIR" "$LOCKDIR"
cd "$PAN" || exit 2
echo "MONITOR_START $(date '+%F %T %Z') node=$NODE_LABEL host=$(hostname) list=$LIST_FILE" >> "$LOG"

active_for_base() {
  base="$1"
  ps -u $USER -o cmd 2>/dev/null \
    | grep -F "$base.mod.panEDTA" \
    | grep -E 'RepeatMasker|rmblastn|ProcessRepeats' \
    | grep -v grep \
    | wc -l
}

all_finished=0
while [ "$all_finished" -eq 0 ]; do
  all_finished=1
  while IFS= read -r genome; do
    [ -z "$genome" ] && continue
    base="$(basename "$genome")"
    out="$PAN/$base.mod.panEDTA.out"
    lock="$LOCKDIR/$base.lock"
    done_marker="$LOGDIR/$base.done"
    failed_marker="$LOGDIR/$base.failed"
    active="$(active_for_base "$base")"

    if [ -s "$done_marker" ] || [ -s "$failed_marker" ]; then
      continue
    fi

    all_finished=0
    echo "$(date '+%F %T') base=$base active=$active out_size=$(stat -c %s "$out" 2>/dev/null || echo 0)" >> "$LOG"

    if [ "$active" -eq 0 ] && [ -s "$out" ]; then
      path/to/home/anaconda3/envs/EDTA_env/bin/perl -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' "$out" 2>> "$LOG" || true
      touch "$done_marker"
      rmdir "$lock" 2>/dev/null || true
      echo "$(date '+%F %T') DONE_CLEAN $base" >> "$LOG"
    elif [ "$active" -eq 0 ] && [ ! -s "$out" ]; then
      touch "$failed_marker"
      rmdir "$lock" 2>/dev/null || true
      echo "$(date '+%F %T') FAILED_NO_OUT $base" >> "$LOG"
    fi
  done < "$LIST_FILE"
  [ "$all_finished" -eq 1 ] && break
  sleep 300
done

echo "MONITOR_END $(date '+%F %T %Z') node=$NODE_LABEL" >> "$LOG"
