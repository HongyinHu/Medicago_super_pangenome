#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
LIST=$BASE/input/valid_phenotyped_bams.for_delly.tsv
export PATH=path/to/home/anaconda3/envs/smoove/bin:$PATH
TASK_ID=${SLURM_ARRAY_TASK_ID:-${1:-}}
if [[ -z "$TASK_ID" ]]; then echo "missing task id" >&2; exit 2; fi
line=$(sed -n "${TASK_ID}p" "$LIST")
if [[ -z "$line" ]]; then echo "no input line for task $TASK_ID" >&2; exit 2; fi
sid=$(echo "$line" | awk -F'\t' '{print $1}')
phe=$(echo "$line" | awk -F'\t' '{print $2}')
bam=$(echo "$line" | awk -F'\t' '{print $3}')
outdir=$BASE/results/smoove_per_sample/$sid
outvcf=$outdir/$sid-smoove.genotyped.vcf.gz
done=$outdir/$sid.smoove.done
if [[ -s "$done" && -s "$outvcf" && -s "$outvcf.tbi" ]]; then
  echo "[$(date '+%F %T')] skip completed SID=$sid"
  exit 0
fi
rm -rf "$outdir"
mkdir -p "$outdir"
echo "[$(date '+%F %T')] smoove start SID=$sid phenotype=$phe host=$(hostname)"
smoove call --outdir "$outdir" --name "$sid" --fasta "$REF" --genotype -p ${SMOOVE_THREADS:-4} "$bam"
printf 'sample_id\t%s\nphenotype\t%s\nbam\t%s\nvcf\t%s\ntime\t%s\n' "$sid" "$phe" "$bam" "$outvcf" "$(date '+%F %T %Z')" > "$done"
echo "[$(date '+%F %T')] smoove done SID=$sid"
