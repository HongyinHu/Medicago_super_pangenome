#!/usr/bin/env bash
set -euo pipefail
RUN=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
DELLY=path/to/home/anaconda3/envs/delly/bin/delly
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
TASK=${SLURM_ARRAY_TASK_ID:-${1:-}}
[ -n "$TASK" ] || { echo "Need task id" >&2; exit 2; }
LINE=$(sed -n "${TASK}p" "$RUN/input/valid_phenotyped_bams.for_delly.tsv")
[ -n "$LINE" ] || { echo "No line for task $TASK" >&2; exit 2; }
SID=$(echo "$LINE" | cut -f1)
PHEN=$(echo "$LINE" | cut -f2)
BAM=$(echo "$LINE" | cut -f3)
OUTDIR="$RUN/results/delly_DEL_per_sample/$SID"
mkdir -p "$OUTDIR"
OUTBCF="$OUTDIR/$SID.delly.DEL.bcf"
OUTVCF="$OUTDIR/$SID.delly.DEL.vcf.gz"
DONE="$OUTDIR/$SID.delly.DEL.done"
if [ -s "$OUTVCF" ] && [ -s "$OUTVCF.tbi" ] && [ -s "$DONE" ]; then
  echo "[$(date '+%F %T')] skip done $SID"
  exit 0
fi
echo "[$(date '+%F %T')] start SID=$SID phenotype=$PHEN bam=$BAM host=$(hostname)"
"$DELLY" call -t DEL -g "$REF" -o "$OUTBCF" "$BAM"
"$BCFTOOLS" view -Oz -o "$OUTVCF" "$OUTBCF"
"$BCFTOOLS" index -t "$OUTVCF"
"$BCFTOOLS" view -r Chr4:88970000-88995000 "$OUTVCF" > "$OUTDIR/$SID.delly.DEL.Chr23997_window.vcf" || true
{
  echo "sample_id=$SID"
  echo "phenotype=$PHEN"
  echo "bam=$BAM"
  echo "finished=$(date '+%F %T')"
  echo "records_total=$($BCFTOOLS view -H "$OUTVCF" | wc -l)"
  echo "records_Chr23997_window=$(grep -vc '^#' "$OUTDIR/$SID.delly.DEL.Chr23997_window.vcf" 2>/dev/null || echo 0)"
} > "$DONE"
echo "[$(date '+%F %T')] done $SID"
