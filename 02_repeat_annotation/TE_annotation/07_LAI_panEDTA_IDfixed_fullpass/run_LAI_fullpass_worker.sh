#!/usr/bin/env bash
set -uo pipefail

BASE=path/to/project/N_1.EDTA_single
OUTDIR=$BASE/07_LAI_panEDTA_IDfixed_fullpass
LOGROOT=$OUTDIR/logs
QUEUE=${QUEUE:-$OUTDIR/genome.queue}
THREADS=${THREADS:-8}

export PATH=path/to/home/anaconda3/envs/EDTA_env/bin:path/to/home/anaconda3/envs/EDTA_env/share/EDTA:path/to/home/anaconda3/envs/EDTA_env/share/LTR_retriever:$PATH
export PYTHONNOUSERSITE=1

mkdir -p "$LOGROOT"/{locks,done,not_applicable,failed,status}

status_line() {
  local sample="$1" lai_file="$2" value="$3" intact="$4" panout="$5" status="$6" note="$7"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$lai_file" "$value" "$intact" "$panout" "$status" "$note" > "$LOGROOT/status/$sample.tsv"
}

process_one() {
  local genome="$1"
  local sample="${genome%.fa}"
  local work="$OUTDIR/$sample"
  local lock="$LOGROOT/locks/$sample.lockdir"
  local final="$work/$genome.mod.panEDTA.IDfixed.out.LAI"
  local stamp stdout stderr rc

  if [ -s "$LOGROOT/done/$sample.done" ] || [ -s "$LOGROOT/not_applicable/$sample.done" ]; then
    return 0
  fi
  if [ -s "$final" ]; then
    status_line "$sample" "$final" NA "$(readlink -f "$work/$genome.mod.pass.list" 2>/dev/null || true)" "$(readlink -f "$work/$genome.mod.panEDTA.IDfixed.out" 2>/dev/null || true)" done existing_final
    touch "$LOGROOT/done/$sample.done"
    return 0
  fi
  if ! mkdir "$lock" 2>/dev/null; then
    return 0
  fi

  stamp=$(date +%Y%m%d_%H%M%S)
  stdout="$LOGROOT/$sample.stdout.$stamp.log"
  stderr="$LOGROOT/$sample.stderr.$stamp.log"

  if [ ! -s "$work/$genome.mod" ] || [ ! -s "$work/$genome.mod.pass.list" ] || [ ! -s "$work/$genome.mod.panEDTA.IDfixed.out" ]; then
    status_line "$sample" "$final" NA "$(readlink -f "$work/$genome.mod.pass.list" 2>/dev/null || true)" "$(readlink -f "$work/$genome.mod.panEDTA.IDfixed.out" 2>/dev/null || true)" failed missing_input
    touch "$LOGROOT/failed/$sample.failed"
    rmdir "$lock" 2>/dev/null || true
    return 2
  fi

  if compgen -G "$work/$genome.mod.panEDTA.IDfixed.out.LAI*" >/dev/null; then
    local backup="$work/recover_backup_$stamp"
    mkdir -p "$backup"
    find "$work" -maxdepth 1 -type f -name "$genome.mod.panEDTA.IDfixed.out.LAI*" -exec mv -t "$backup" {} + 2>/dev/null || true
  fi

  echo "[$(date)] RUN $sample host=$(hostname) threads=$THREADS" | tee -a "$LOGROOT/driver.log"
  (
    cd "$work" || exit 2
    LAI -genome "$genome.mod" -intact "$genome.mod.pass.list" -all "$genome.mod.panEDTA.IDfixed.out" -t "$THREADS" > "$stdout" 2> "$stderr"
  )
  rc=$?

  if [ -s "$final" ]; then
    status_line "$sample" "$final" NA "$(readlink -f "$work/$genome.mod.pass.list" 2>/dev/null || true)" "$(readlink -f "$work/$genome.mod.panEDTA.IDfixed.out" 2>/dev/null || true)" done "rc=$rc"
    touch "$LOGROOT/done/$sample.done"
    rmdir "$lock" 2>/dev/null || true
    return 0
  fi
  if grep -q 'LAI is not applicable' "$stdout" 2>/dev/null; then
    status_line "$sample" "$final" NA "$(readlink -f "$work/$genome.mod.pass.list" 2>/dev/null || true)" "$(readlink -f "$work/$genome.mod.panEDTA.IDfixed.out" 2>/dev/null || true)" not_applicable low_intact
    touch "$LOGROOT/not_applicable/$sample.done"
    rmdir "$lock" 2>/dev/null || true
    return 0
  fi

  status_line "$sample" "$final" NA "$(readlink -f "$work/$genome.mod.pass.list" 2>/dev/null || true)" "$(readlink -f "$work/$genome.mod.panEDTA.IDfixed.out" 2>/dev/null || true)" failed "rc=$rc"
  touch "$LOGROOT/failed/$sample.failed"
  rmdir "$lock" 2>/dev/null || true
  return "$rc"
}

rc=0
while IFS= read -r genome; do
  [ -n "$genome" ] || continue
  process_one "$genome" || rc=1
done < "$QUEUE"
exit "$rc"
