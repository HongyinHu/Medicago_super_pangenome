#!/usr/bin/env bash
set -euo pipefail
OLD=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706
MED=$OLD/medicago144_3caller_20260709
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
CONFIG=path/to/home/anaconda3/envs/manta/share/manta-1.6.0-2/bin/configManta.py
TASK=${SLURM_ARRAY_TASK_ID:-${1:-}}
[ -n "$TASK" ] || { echo "Need task id" >&2; exit 2; }
LINE=$(awk -F'\t' -v n="$TASK" 'NR==n+1{print}' "$MED/input/medicago144_sample_bams.tsv")
[ -n "$LINE" ] || { echo "No line for task $TASK" >&2; exit 2; }
SID=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE" | cut -f6)
RUNDIR="$MED/results/manta_per_sample/$SID"
LOCKDIR="$MED/locks/manta.$SID.lock"
if ! mkdir "$LOCKDIR" 2>/dev/null; then
  echo "[$(date '+%F %T')] manta lock exists, skip SID=$SID"
  exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT
printf 'host=%s
pid=%s
time=%s
' "$(hostname)" "$$" "$(date '+%F %T %Z')" > "$LOCKDIR/info.txt"
OUTVCF="$RUNDIR/results/variants/diploidSV.vcf.gz"
DONE="$RUNDIR/$SID.manta.done"
FAIL="$RUNDIR/$SID.manta.failed"
if [ -s "$DONE" ] && [ -s "$OUTVCF" ] && [ -s "$OUTVCF.tbi" ]; then echo "skip manta $SID"; exit 0; fi
OLDVCF="$OLD/results/manta_per_sample/$SID/results/variants/diploidSV.vcf.gz"
OLDDONE="$OLD/results/manta_per_sample/$SID/$SID.manta.done"
if [ -s "$OLDVCF" ] && [ -s "$OLDVCF.tbi" ] && [ -s "$OLDDONE" ]; then
  mkdir -p "$RUNDIR/results/variants"
  ln -sf "$OLDVCF" "$OUTVCF"
  ln -sf "$OLDVCF.tbi" "$OUTVCF.tbi"
  printf 'sample_id\t%s\nreused_from\t%s\ntime\t%s\n' "$SID" "$OLDVCF" "$(date '+%F %T %Z')" > "$DONE"
  echo "reuse old manta $SID"
  exit 0
fi
if grep -q "^$SID\b" "$OLD/summary/manta_failed_samples.tsv" 2>/dev/null; then
  mkdir -p "$RUNDIR"
  grep "^$SID\b" "$OLD/summary/manta_failed_samples.tsv" > "$FAIL"
  echo "known old manta failure $SID; mark unavailable"
  exit 0
fi
rm -rf "$RUNDIR"
mkdir -p "$RUNDIR"
echo "[$(date '+%F %T')] Manta start SID=$SID host=$(hostname)"
set +e
"$CONFIG" --bam "$BAM" --referenceFasta "$REF" --runDir "$RUNDIR" > "$RUNDIR/configManta.stdout" 2> "$RUNDIR/configManta.stderr"
rc1=$?
if [ "$rc1" -eq 0 ]; then
  "$RUNDIR/runWorkflow.py" -m local -j ${MANTA_THREADS:-4} > "$RUNDIR/runWorkflow.stdout" 2> "$RUNDIR/runWorkflow.stderr"
  rc2=$?
else
  rc2=999
fi
set -e
if [ "$rc1" -ne 0 ] || [ "$rc2" -ne 0 ] || [ ! -s "$OUTVCF" ]; then
  { echo -e "sample_id\t$SID"; echo -e "bam\t$BAM"; echo -e "config_rc\t$rc1"; echo -e "workflow_rc\t$rc2"; echo -e "time\t$(date '+%F %T %Z')"; tail -80 "$RUNDIR/runWorkflow.stderr" 2>/dev/null || true; } > "$FAIL"
  echo "Manta failed SID=$SID rc=$rc1/$rc2; marked unavailable"
  exit 0
fi
printf 'sample_id\t%s\nbam\t%s\nvcf\t%s\ntime\t%s\n' "$SID" "$BAM" "$OUTVCF" "$(date '+%F %T %Z')" > "$DONE"
echo "[$(date '+%F %T')] Manta done SID=$SID"
