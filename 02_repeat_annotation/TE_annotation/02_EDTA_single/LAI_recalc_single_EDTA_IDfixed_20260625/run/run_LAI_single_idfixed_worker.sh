#!/usr/bin/env bash
set -uo pipefail
BASE=path/to/project/N_1.EDTA_single
OUT=$BASE/02_EDTA_single/LAI_recalc_single_EDTA_IDfixed_20260625
LOGROOT=$OUT/logs
QUEUE=${QUEUE:-$OUT/sample.queue}
THREADS=${THREADS:-30}
NODE=${NODE:-$(hostname -s)}
export PATH=path/to/home/anaconda3/envs/EDTA_env/bin:path/to/home/anaconda3/envs/EDTA_env/share/EDTA:path/to/home/anaconda3/envs/EDTA_env/share/LTR_retriever:$PATH
export PYTHONNOUSERSITE=1
mkdir -p "$LOGROOT"/{locks,done,not_applicable,failed,status} "$LOGROOT/$NODE"
status_line() {
  local sample="$1" lai_file="$2" value="$3" status="$4" note="$5"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$lai_file" "$value" "$status" "$note" "$(date '+%F %T')" > "$LOGROOT/status/$sample.tsv"
}
get_field() { awk -F '\t' -v s="$1" -v col="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $1==s{print $h[col]; exit}' "$OUT/prepare_manifest_single_idfixed.tsv"; }
process_one() {
  local sample="$1" work final lock stdout stderr rc value
  work=$(get_field "$sample" workdir)
  final="$work/$sample.singleEDTA.IDfixed.out.LAI"
  lock="$LOGROOT/locks/$sample.lockdir"
  [ -s "$LOGROOT/done/$sample.done" ] || [ -s "$LOGROOT/not_applicable/$sample.done" ] && return 0
  [ -s "$LOGROOT/failed/$sample.failed" ] && return 0
  if [ -s "$final" ]; then
    value=$(awk '$1=="whole_genome"{print $7; exit}' "$final" 2>/dev/null || echo NA)
    status_line "$sample" "$final" "${value:-NA}" done existing_final
    touch "$LOGROOT/done/$sample.done"
    return 0
  fi
  if ! mkdir "$lock" 2>/dev/null; then return 0; fi
  echo -e "$NODE\t$(hostname)\t$(date '+%F %T')" > "$lock/info.tsv"
  stdout="$LOGROOT/$NODE/$sample.stdout.$(date +%Y%m%d_%H%M%S).log"
  stderr="$LOGROOT/$NODE/$sample.stderr.$(date +%Y%m%d_%H%M%S).log"
  echo "[$(date)] RUN $sample host=$(hostname) node=$NODE threads=$THREADS" >> "$LOGROOT/driver.log"
  (
    cd "$work" || exit 2
    rm -f "$sample.singleEDTA.IDfixed.out.LAI" "$sample.singleEDTA.IDfixed.out.LAI".* 2>/dev/null || true
    LAI -genome "$sample.singleEDTA.genome.fa" -intact "$sample.singleEDTA.pass.list" -all "$sample.singleEDTA.IDfixed.out" -t "$THREADS" > "$stdout" 2> "$stderr"
  )
  rc=$?
  if [ -s "$final" ]; then
    value=$(awk '$1=="whole_genome"{print $7; exit}' "$final" 2>/dev/null || echo NA)
    status_line "$sample" "$final" "${value:-NA}" done "rc=$rc"
    touch "$LOGROOT/done/$sample.done"
    rm -f "$LOGROOT/failed/$sample.failed"
    rmdir "$lock" 2>/dev/null || true
    return 0
  fi
  if grep -q 'LAI is not applicable' "$stdout" "$stderr" 2>/dev/null; then
    status_line "$sample" "$final" NA not_applicable low_content_or_no_ltr
    touch "$LOGROOT/not_applicable/$sample.done"
    rmdir "$lock" 2>/dev/null || true
    return 0
  fi
  status_line "$sample" "$final" NA failed "rc=$rc stdout=$stdout stderr=$stderr"
  touch "$LOGROOT/failed/$sample.failed"
  rmdir "$lock" 2>/dev/null || true
  return 0
}
rc=0
while IFS= read -r sample; do
  [ -n "$sample" ] || continue
  process_one "$sample" || rc=1
done < "$QUEUE"
exit "$rc"
