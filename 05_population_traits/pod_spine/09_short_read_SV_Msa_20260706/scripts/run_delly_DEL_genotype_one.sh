#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
DELLY=path/to/home/anaconda3/envs/delly/bin/delly
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
SITES=$BASE/results/delly_DEL_cohort/phenotyped99.DEL.sites.bcf
LIST=$BASE/input/valid_phenotyped_bams.for_delly.tsv
TASK_ID=${SLURM_ARRAY_TASK_ID:-${1:-}}
if [[ -z "$TASK_ID" ]]; then echo "missing task id" >&2; exit 2; fi
line=$(sed -n "${TASK_ID}p" "$LIST")
if [[ -z "$line" ]]; then echo "no input line for task $TASK_ID" >&2; exit 2; fi
sid=$(echo "$line" | awk -F'\t' '{print $1}')
phe=$(echo "$line" | awk -F'\t' '{print $2}')
bam=$(echo "$line" | awk -F'\t' '{print $3}')
outdir=$BASE/results/delly_DEL_genotype_per_sample/$sid
mkdir -p "$outdir"
outbcf=$outdir/$sid.delly.DEL.geno.bcf
done=$outdir/$sid.delly.DEL.geno.done
if [[ -s "$done" && -s "$outbcf" && -s "$outbcf.csi" ]]; then
  echo "[$(date '+%F %T')] skip completed SID=$sid"
  exit 0
fi
rm -f "$outbcf" "$outbcf.csi" "$done"
echo "[$(date '+%F %T')] genotype start SID=$sid phenotype=$phe host=$(hostname)"
"$DELLY" call -g "$REF" -v "$SITES" -o "$outbcf" "$bam"
"$BCFTOOLS" index -f "$outbcf"
printf 'sample_id\t%s\nphenotype\t%s\nbam\t%s\ntime\t%s\n' "$sid" "$phe" "$bam" "$(date '+%F %T %Z')" > "$done"
echo "[$(date '+%F %T')] genotype done SID=$sid"
