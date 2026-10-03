#!/usr/bin/env bash
set -euo pipefail
OLD=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706
MED=$OLD/medicago144_3caller_20260709
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
export PATH=path/to/home/anaconda3/envs/smoove/bin:$PATH
TASK=${SLURM_ARRAY_TASK_ID:-${1:-}}
[ -n "$TASK" ] || { echo "Need task id" >&2; exit 2; }
LINE=$(awk -F'\t' -v n="$TASK" 'NR==n+1{print}' "$MED/input/medicago144_sample_bams.tsv")
[ -n "$LINE" ] || { echo "No line for task $TASK" >&2; exit 2; }
SID=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE" | cut -f6)
OUTDIR="$MED/results/smoove_per_sample/$SID"
LOCKDIR="$MED/locks/smoove.$SID.lock"
if ! mkdir "$LOCKDIR" 2>/dev/null; then
  echo "[$(date '+%F %T')] smoove lock exists, skip SID=$SID"
  exit 0
fi
trap 'rm -rf "$LOCKDIR"' EXIT
printf 'host=%s
pid=%s
time=%s
' "$(hostname)" "$$" "$(date '+%F %T %Z')" > "$LOCKDIR/info.txt"
OUTVCF="$OUTDIR/$SID-smoove.genotyped.vcf.gz"
DONE="$OUTDIR/$SID.smoove.done"
FAIL="$OUTDIR/$SID.smoove.failed"
if [ -s "$DONE" ] && [ -s "$OUTVCF" ] && [ -s "$OUTVCF.tbi" ]; then echo "skip smoove $SID"; exit 0; fi
OLDVCF="$OLD/results/smoove_per_sample/$SID/$SID-smoove.genotyped.vcf.gz"
OLDDONE="$OLD/results/smoove_per_sample/$SID/$SID.smoove.done"
if [ -s "$OLDVCF" ] && [ -s "$OLDVCF.tbi" ] && [ -s "$OLDDONE" ]; then
  mkdir -p "$OUTDIR"
  ln -sf "$OLDVCF" "$OUTVCF"
  ln -sf "$OLDVCF.tbi" "$OUTVCF.tbi"
  printf 'sample_id\t%s\nreused_from\t%s\ntime\t%s\n' "$SID" "$OLDVCF" "$(date '+%F %T %Z')" > "$DONE"
  echo "reuse old smoove $SID"
  exit 0
fi
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"
echo "[$(date '+%F %T')] smoove start SID=$SID host=$(hostname)"
set +e
smoove call --outdir "$OUTDIR" --name "$SID" --fasta "$REF" --genotype -p ${SMOOVE_THREADS:-4} "$BAM" > "$OUTDIR/smoove.stdout" 2> "$OUTDIR/smoove.stderr"
rc=$?
set -e
if [ "$rc" -ne 0 ] || [ ! -s "$OUTVCF" ]; then
  { echo -e "sample_id\t$SID"; echo -e "bam\t$BAM"; echo -e "rc\t$rc"; echo -e "time\t$(date '+%F %T %Z')"; tail -80 "$OUTDIR/smoove.stderr" 2>/dev/null || true; } > "$FAIL"
  echo "smoove failed SID=$SID rc=$rc; marked unavailable"
  exit 0
fi
printf 'sample_id\t%s\nbam\t%s\nvcf\t%s\ntime\t%s\n' "$SID" "$BAM" "$OUTVCF" "$(date '+%F %T %Z')" > "$DONE"
echo "[$(date '+%F %T')] smoove done SID=$SID"
