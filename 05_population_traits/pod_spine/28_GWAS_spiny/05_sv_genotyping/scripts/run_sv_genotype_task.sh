#!/usr/bin/env bash
set -euo pipefail

RUN=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$RUN/05_sv_genotyping"
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
MANIFEST="$STAGE/summary/sv_analysis_manifest.tsv"
SITES="$STAGE/candidate_sites/medicago143.sites.filtered.nsamp2.smoove.vcf.gz"
SITES_INDEXED="$STAGE/candidate_sites/medicago143.sites.filtered.nsamp2.vcf.gz"
SMOOVE=path/to/home/anaconda3/envs/smoove/bin/smoove
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}

export PATH=path/to/home/anaconda3/envs/smoove/bin:path/to/home/anaconda3/envs/panpop/bin:/usr/bin:/bin
line=$(awk -F '\t' -v n=$((TASK_ID + 1)) 'NR==n {print; exit}' "$MANIFEST")
[[ -n "$line" ]] || { echo "No manifest row for task $TASK_ID" >&2; exit 1; }
IFS=$'\t' read -r sid bam discovery_vcf filtered_vcf <<< "$line"
OUTDIR="$STAGE/joint_genotype_per_sample/$sid"
OUT="$OUTDIR/$sid-smoove.genotyped.vcf.gz"
DONE="$OUTDIR/$sid.done"
mkdir -p "$OUTDIR" "$STAGE/logs/genotype"
EXPECTED_RECORDS=$("$BCFTOOLS" index -n "$SITES_INDEXED")

if [[ -s "$OUT" && -s "$OUT.csi" ]] \
    && "$BCFTOOLS" view -h "$OUT" >/dev/null 2>&1 \
    && [[ "$("$BCFTOOLS" index -n "$OUT")" == "$EXPECTED_RECORDS" ]] \
    && [[ "$("$BCFTOOLS" query -l "$OUT")" == "$sid" ]]; then
    if [[ ! -s "$DONE" ]]; then
        printf 'completed\t%s\nhost\t%s\nsample\t%s\nrecords\t%s\nmode\trecovered_valid_genotype_no_duphold\n' \
            "$(date '+%F %T %Z')" "$(hostname)" "$sid" "$EXPECTED_RECORDS" > "$DONE"
    fi
    echo "$sid already complete"
    exit 0
fi
rm -f "$DONE"
"$SMOOVE" genotype --name "$sid" --outdir "$OUTDIR" --fasta "$REF" \
    --removepr --processes "${SLURM_CPUS_PER_TASK:-2}" --vcf "$SITES" "$bam"
[[ -s "$OUT" ]]
[[ -s "$OUT.csi" ]]
"$BCFTOOLS" view -h "$OUT" >/dev/null
[[ "$("$BCFTOOLS" index -n "$OUT")" == "$EXPECTED_RECORDS" ]]
[[ "$("$BCFTOOLS" query -l "$OUT")" == "$sid" ]]
printf 'completed\t%s\nhost\t%s\nsample\t%s\nrecords\t%s\nmode\tgenotype_no_duphold\n' \
    "$(date '+%F %T %Z')" "$(hostname)" "$sid" "$EXPECTED_RECORDS" > "$DONE"
