#!/usr/bin/env bash
set -euo pipefail
MED=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706/medicago144_3caller_20260709
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
DELLY=path/to/home/anaconda3/envs/delly/bin/delly
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
TASK=${SLURM_ARRAY_TASK_ID:-${1:-}}
[ -n "$TASK" ] || { echo "Need task id" >&2; exit 2; }
LINE=$(awk -F'\t' -v n="$TASK" 'NR==n+1{print}' "$MED/input/medicago144_sample_bams.tsv")
[ -n "$LINE" ] || { echo "No line for task $TASK" >&2; exit 2; }
SID=$(echo "$LINE" | cut -f1)
BAM=$(echo "$LINE" | cut -f6)
OUTDIR="$MED/results/delly_per_sample/$SID"
mkdir -p "$OUTDIR"
OUTBCF="$OUTDIR/$SID.delly.all.bcf"
OUTVCF="$OUTDIR/$SID.delly.all.vcf.gz"
DONE="$OUTDIR/$SID.delly.all.done"
if [ -s "$DONE" ] && [ -s "$OUTVCF" ] && [ -s "$OUTVCF.tbi" ]; then echo "skip delly $SID"; exit 0; fi
rm -f "$OUTBCF" "$OUTBCF.csi" "$OUTVCF" "$OUTVCF.tbi" "$DONE"
echo "[$(date '+%F %T')] DELLY start SID=$SID host=$(hostname)"
"$DELLY" call -g "$REF" -o "$OUTBCF" "$BAM"
"$BCFTOOLS" index -f "$OUTBCF"
"$BCFTOOLS" view -Oz -o "$OUTVCF" "$OUTBCF"
"$BCFTOOLS" index -t -f "$OUTVCF"
printf 'sample_id\t%s\nbam\t%s\nvcf\t%s\ntime\t%s\n' "$SID" "$BAM" "$OUTVCF" "$(date '+%F %T %Z')" > "$DONE"
echo "[$(date '+%F %T')] DELLY done SID=$SID"
