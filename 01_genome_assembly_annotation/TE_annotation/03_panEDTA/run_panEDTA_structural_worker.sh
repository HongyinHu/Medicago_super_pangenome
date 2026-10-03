#!/usr/bin/env bash
set -uo pipefail

PAN="path/to/project/N_1.EDTA_single/03_panEDTA"
LIST="${1:?list file required}"
THREADS="${THREADS:-12}"
NODE_NAME="${NODE_NAME:-$(hostname)}"
LOGDIR="$PAN/logs/structural_resume/$NODE_NAME"
mkdir -p "$LOGDIR" "$PAN/logs/structural_resume/reconciled"
cd "$PAN" || exit 2

source ~/.bashrc >> "$LOGDIR/driver.log" 2>&1 || true
conda activate EDTA_env >> "$LOGDIR/driver.log" 2>&1 || {
  echo "ERROR $(date '+%F %T') conda activate EDTA_env failed" >> "$LOGDIR/driver.log"
  exit 2
}
export PYTHONNOUSERSITE=1
EDTA_PL="$(command -v EDTA.pl || true)"
if [ ! -x "$EDTA_PL" ]; then
  echo "ERROR $(date '+%F %T') EDTA.pl not found in EDTA_env" >> "$LOGDIR/driver.log"
  exit 2
fi

{
  echo "START $(date '+%F %T %Z') node=$NODE_NAME host=$(hostname) threads=$THREADS list=$LIST"
  echo "EDTA=$EDTA_PL"
} >> "$LOGDIR/driver.log"

while IFS= read -r line || [ -n "$line" ]; do
  genome_path="$(echo "$line" | awk '{print $1}')"
  cds_ind="$(echo "$line" | awk '{print $2}')"
  genome="$(basename "$genome_path" 2>/dev/null)"
  [ -z "$genome" ] && continue

  if [ -s "logs/structural_resume/reconciled/$genome.done" ]; then
    echo "SKIP $(date '+%F %T') $genome reconciled_done" >> "$LOGDIR/driver.log"
    continue
  fi
  if [ ! -s "$genome.mod.panEDTA.out" ]; then
    echo "FAILED $(date '+%F %T') $genome missing $genome.mod.panEDTA.out" >> "$LOGDIR/driver.log"
    echo "missing_rmout" > "$LOGDIR/$genome.failed"
    continue
  fi

  perl -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' "$genome.mod.panEDTA.out"

  cds_args=()
  if [ -n "$cds_ind" ]; then
    cds_base="$(basename "$cds_ind" 2>/dev/null)"
    if [ -s "$cds_base" ]; then
      cds_args=(--cds "$cds_base")
    elif [ -s "$cds_ind" ]; then
      cds_args=(--cds "$cds_ind")
    else
      echo "WARN $(date '+%F %T') $genome cds_not_found=$cds_ind" >> "$LOGDIR/driver.log"
    fi
  fi

  glog="$LOGDIR/$genome.EDTA_final.$(date '+%Y%m%d_%H%M%S').log"
  echo "RUN $(date '+%F %T') $genome log=$glog" >> "$LOGDIR/driver.log"
  perl "$EDTA_PL" --genome "$genome" -t "$THREADS" --step final --anno 1 --curatedlib genome.list.panEDTA.TElib.fa "${cds_args[@]}" --rmout "$genome.mod.panEDTA.out" > "$glog" 2>&1
  rc=$?
  if [ "$rc" -eq 0 ] && [ -s "$genome.mod.EDTA.TEanno.sum" ] && grep -q 'Evaluation of TE annotation finished' "$glog"; then
    echo "DONE $(date '+%F %T') $genome sum_size=$(stat -c %s "$genome.mod.EDTA.TEanno.sum")" >> "$LOGDIR/driver.log"
    rm -f logs/structural_resume/*/$genome.failed 2>/dev/null || true
    echo "done $(date '+%F %T') sum=$genome.mod.EDTA.TEanno.sum log=$glog" > "logs/structural_resume/reconciled/$genome.done"
  else
    echo "FAILED $(date '+%F %T') $genome rc=$rc expected_success_log_or_sum_missing" >> "$LOGDIR/driver.log"
    echo "rc=$rc" > "$LOGDIR/$genome.failed"
  fi
done < "$LIST"

echo "END $(date '+%F %T %Z') node=$NODE_NAME" >> "$LOGDIR/driver.log"