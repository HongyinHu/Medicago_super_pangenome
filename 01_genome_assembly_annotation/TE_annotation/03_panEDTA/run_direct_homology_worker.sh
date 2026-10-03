#!/usr/bin/env bash
set -uo pipefail

PAN="path/to/project/N_1.EDTA_single/03_panEDTA"
LIST_FILE="${1:?list file required}"
NODE_LABEL="${2:?node label required}"
THREADS_PER_JOB="${3:-16}"
MAX_JOBS="${4:-1}"
MAX_LOAD="${5:-100}"

LOGDIR="$PAN/logs/direct_homology/$NODE_LABEL"
LOCKDIR="$PAN/logs/homology_locks"
LIB="$PAN/genome.list.panEDTA.TElib.fa"

mkdir -p "$LOGDIR" "$LOCKDIR"
cd "$PAN" || exit 2

EDTA_ENV_BIN="path/to/home/anaconda3/envs/EDTA_env/bin"
export PATH="$EDTA_ENV_BIN:$PATH"
export PYTHONNOUSERSITE=1
if ! command -v RepeatMasker >/dev/null 2>&1; then
  echo "ERROR: RepeatMasker not found after adding $EDTA_ENV_BIN" >&2
  exit 3
fi

if [ ! -s "$LIB" ]; then
  echo "ERROR: missing panEDTA library $LIB" >&2
  exit 4
fi

echo "START $(date '+%F %T %Z') node=$NODE_LABEL host=$(hostname) threads=$THREADS_PER_JOB max_jobs=$MAX_JOBS max_load=$MAX_LOAD list=$LIST_FILE"
echo "RepeatMasker=$(command -v RepeatMasker || echo NA)"

current_load_int() {
  awk '{split($1,a,"."); print a[1]}' /proc/loadavg 2>/dev/null || echo 0
}

wait_for_load() {
  while true; do
    load_now="$(current_load_int)"
    if [ "${load_now:-0}" -lt "$MAX_LOAD" ]; then
      return 0
    fi
    echo "$(date '+%F %T') load=$load_now >= $MAX_LOAD; waiting"
    sleep 300
  done
}

run_one() {
  genome="$1"
  base="$(basename "$genome")"
  mod="$PAN/$base.mod"
  target="$PAN/$base.mod.panEDTA"
  out="$PAN/$base.mod.panEDTA.out"
  lock="$LOCKDIR/$base.lock"
  marker_done="$LOGDIR/$base.done"
  marker_failed="$LOGDIR/$base.failed"
  log="$LOGDIR/$base.RepeatMasker.$(date +%Y%m%d_%H%M%S).log"

  if [ -s "$out" ]; then
    echo "$(date '+%F %T') SKIP existing $out"
    touch "$marker_done"
    return 0
  fi

  if ! mkdir "$lock" 2>/dev/null; then
    echo "$(date '+%F %T') SKIP locked $base"
    return 0
  fi

  rm -f "$marker_failed"
  {
    echo "BEGIN $(date '+%F %T %Z') node=$NODE_LABEL host=$(hostname) genome=$genome base=$base"
    wait_for_load
    if [ ! -e "$mod" ] && [ -e "$genome.mod" ]; then
      ln -sf "$genome.mod" "$mod"
    fi
    if [ ! -e "$mod" ]; then
      echo "ERROR: missing mod fasta $mod and $genome.mod"
      exit 10
    fi
    ln -sf "$base.mod" "$target"
    RepeatMasker -e ncbi -pa "$THREADS_PER_JOB" -q -div 40 -lib "$LIB" -cutoff 225 -gff "$target"
    rc=$?
    if [ "$rc" -eq 0 ] && [ -s "$out" ]; then
      perl -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' "$out"
      echo "DONE $(date '+%F %T %Z') out=$out size=$(stat -c %s "$out" 2>/dev/null || echo NA)"
      exit 0
    fi
    echo "ERROR: RepeatMasker rc=$rc out_missing_or_empty=$out"
    exit "${rc:-1}"
  } > "$log" 2>&1
  rc=$?

  rmdir "$lock" 2>/dev/null || true
  if [ "$rc" -eq 0 ]; then
    touch "$marker_done"
    return 0
  fi
  touch "$marker_failed"
  echo "$(date '+%F %T') FAILED $base rc=$rc log=$log" >&2
  return "$rc"
}

while IFS= read -r genome; do
  [ -z "$genome" ] && continue
  while [ "$(jobs -pr | wc -l)" -ge "$MAX_JOBS" ]; do
    sleep 15
  done
  run_one "$genome" &
done < "$LIST_FILE"

while [ "$(jobs -pr | wc -l)" -gt 0 ]; do
  sleep 30
done

fail_count="$(find "$LOGDIR" -maxdepth 1 -name '*.failed' | wc -l)"
done_count="$(find "$LOGDIR" -maxdepth 1 -name '*.done' | wc -l)"
echo "END $(date '+%F %T %Z') node=$NODE_LABEL done_markers=$done_count failed_markers=$fail_count"
[ "$fail_count" -eq 0 ]
