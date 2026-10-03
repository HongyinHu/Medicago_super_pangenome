#!/usr/bin/env bash
set -uo pipefail
PAN="path/to/project/N_1.EDTA_single/03_panEDTA"
QUEUE="${1:?queue required}"
WORKER_ID="${WORKER_ID:-dynamic_fast}"
THREADS="${THREADS:-8}"
NODE_NAME="${NODE_NAME:-$(hostname)}"
EDTA_ENV="path/to/home/anaconda3/envs/EDTA_env"
EDTA_PL="$EDTA_ENV/bin/EDTA.pl"
PERL="$EDTA_ENV/bin/perl"
export PATH="$EDTA_ENV/bin:$PATH"
export PYTHONNOUSERSITE=1
LOGDIR="$PAN/logs/structural_resume/$WORKER_ID"
LOCKDIR="$PAN/logs/structural_resume/locks"
mkdir -p "$LOGDIR" "$LOCKDIR" "$PAN/logs/structural_resume/reconciled"
cd "$PAN" || exit 2
echo "START $(date '+%F %T %Z') worker=$WORKER_ID node=$NODE_NAME host=$(hostname) threads=$THREADS queue=$QUEUE EDTA=$EDTA_PL" >> "$LOGDIR/driver.log"
while true; do
  claimed=""
  while IFS= read -r genome || [ -n "$genome" ]; do
    [ -z "$genome" ] && continue
    [ -s "logs/structural_resume/reconciled/$genome.done" ] && continue
    lock="$LOCKDIR/$genome.lock"
    if mkdir "$lock" 2>/dev/null; then
      echo "$(hostname) $$ $(date '+%F %T %Z')" > "$lock/owner"
      claimed="$genome"
      break
    fi
  done < "$QUEUE"
  [ -z "$claimed" ] && break
  genome="$claimed"
  if [ ! -s "$genome.mod.panEDTA.out" ]; then
    echo "FAILED $(date '+%F %T') $genome missing rmout" >> "$LOGDIR/driver.log"
    echo missing_rmout > "$LOGDIR/$genome.failed"
    rmdir "$LOCKDIR/$genome.lock" 2>/dev/null || rm -rf "$LOCKDIR/$genome.lock"
    continue
  fi
  glog="$LOGDIR/$genome.EDTA_final.$(date '+%Y%m%d_%H%M%S').log"
  echo "RUN $(date '+%F %T') $genome threads=$THREADS log=$glog" >> "$LOGDIR/driver.log"
  "$PERL" -i -nle 's/\s+DNA\s+/\tDNA\/unknown\t/; print $_' "$genome.mod.panEDTA.out"
  "$PERL" "$EDTA_PL" --genome "$genome" -t "$THREADS" --step final --anno 1 --curatedlib genome.list.panEDTA.TElib.fa --rmout "$genome.mod.panEDTA.out" > "$glog" 2>&1
  rc=$?
  if [ "$rc" -eq 0 ] && [ -s "$genome.mod.EDTA.TEanno.sum" ] && grep -q 'Evaluation of TE annotation finished' "$glog"; then
    echo "DONE $(date '+%F %T') $genome" >> "$LOGDIR/driver.log"
    rm -f logs/structural_resume/*/$genome.failed 2>/dev/null || true
    echo "done $(date '+%F %T') sum=$genome.mod.EDTA.TEanno.sum log=$glog" > "logs/structural_resume/reconciled/$genome.done"
  else
    echo "FAILED $(date '+%F %T') $genome rc=$rc" >> "$LOGDIR/driver.log"
    echo "rc=$rc" > "$LOGDIR/$genome.failed"
  fi
  rmdir "$LOCKDIR/$genome.lock" 2>/dev/null || rm -rf "$LOCKDIR/$genome.lock"
done
echo "END $(date '+%F %T %Z') worker=$WORKER_ID" >> "$LOGDIR/driver.log"
