#!/usr/bin/env bash
set -euo pipefail
BASE=${BASE:-path/to/project/N_3.call_SV}
OUT=${OUT:-$BASE/02.assembly_svgap_dualref}
TASKS=${TASKS:-$OUT/01_metadata/svgap_wga_tasks.tsv}
MINIMAP2=${MINIMAP2:-path/to/home/anaconda3/envs/minimap2_env/bin/minimap2}
THREADS=${THREADS:-24}
TASK_ID=${1:?usage: run_svgap_wga_task.sh TASK_ID}
mkdir -p "$OUT/status/wga" "$OUT/logs/tasks" "$OUT/tmp"
line=$(awk -F '\t' -v id="$TASK_ID" 'NR>1 && $1==id{print; found=1; exit} END{if(!found) exit 2}' "$TASKS") || { echo "TASK_NOT_FOUND $TASK_ID" >&2; exit 2; }
IFS=$'\t' read -r task_id ref_run ref_id query_id query_assembly_id source_group role pod_coiling pod_spine independent_trait_sample ref_fasta query_fasta paf <<< "$line"
log="$OUT/logs/tasks/${TASK_ID}.$(hostname).$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$log") 2>&1
fail(){ echo "FAILED task=$TASK_ID time=$(date '+%F %T %Z')"; touch "$OUT/status/wga/${TASK_ID}.failed"; }
trap fail ERR
if [ -s "$paf" ] && [ -e "$OUT/status/wga/${TASK_ID}.done" ]; then
  echo "SKIP_DONE task=$TASK_ID paf=$paf"
  exit 0
fi
rm -f "$OUT/status/wga/${TASK_ID}.failed"
mkdir -p "$(dirname "$paf")"
tmp="${paf}.tmp.$(hostname).$$"
rm -f "$tmp"
echo "START_WGA task=$TASK_ID ref=$ref_id query=$query_id assembly=$query_assembly_id host=$(hostname) threads=$THREADS time=$(date '+%F %T %Z')"
echo "CMD $MINIMAP2 -x asm20 -c --cs=long -t $THREADS $ref_fasta $query_fasta > $paf"
"$MINIMAP2" -x asm20 -c --cs=long -t "$THREADS" "$ref_fasta" "$query_fasta" > "$tmp"
test -s "$tmp"
mv "$tmp" "$paf"
touch "$OUT/status/wga/${TASK_ID}.done"
rm -f "$OUT/status/wga/${TASK_ID}.failed"
echo "DONE_WGA task=$TASK_ID size=$(stat -c %s "$paf") time=$(date '+%F %T %Z')"
