#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/05_sv_genotyping"
STRUCT="$ROOT/03_population_structure"
SMOOVE=path/to/home/anaconda3/envs/smoove/bin/smoove
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools

export PATH=path/to/home/anaconda3/envs/smoove/bin:path/to/home/anaconda3/envs/panpop/bin:/usr/bin:/bin

COHORT="$STAGE/cohort"
TMP="$STAGE/tmp/cohort_qc"
SUMMARY="$STAGE/summary"
SAMPLE_ORDER="$SUMMARY/sv_analysis_samples.order.txt"
VCF_LIST="$SUMMARY/sv_joint_genotype_vcfs.order.txt"
JOINT="$COHORT/medicago143.smoove.square.vcf.gz"
QC="$COHORT/medicago143.smoove.square.gq20.dp5.cr80.maf05.mac5.bnd1.vcf.gz"
EXPECTED_SAMPLES=143
EXPECTED_RECORDS=641339
MIN_CALLED_SAMPLES=115

mkdir -p "$COHORT" "$TMP" "$SUMMARY" "$STAGE/logs/slurm"
test -s "$SUMMARY/sv_joint_genotypes.validated.done"
test -s "$STRUCT/results/plink/analysis.SNP.fam"

awk '{print $2}' "$STRUCT/results/plink/analysis.SNP.fam" > "$SAMPLE_ORDER"
[[ "$(wc -l < "$SAMPLE_ORDER")" -eq "$EXPECTED_SAMPLES" ]]
[[ "$(sort -u "$SAMPLE_ORDER" | wc -l)" -eq "$EXPECTED_SAMPLES" ]]

: > "$VCF_LIST"
while IFS= read -r sid; do
    vcf="$STAGE/joint_genotype_per_sample/$sid/$sid-smoove.genotyped.vcf.gz"
    done_file="$STAGE/joint_genotype_per_sample/$sid/$sid.done"
    test -s "$vcf"
    test -s "$vcf.csi"
    test -s "$done_file"
    [[ "$($BCFTOOLS query -l "$vcf")" == "$sid" ]]
    [[ "$($BCFTOOLS index -n "$vcf")" -eq "$EXPECTED_RECORDS" ]]
    printf '%s\n' "$vcf" >> "$VCF_LIST"
done < "$SAMPLE_ORDER"
[[ "$(wc -l < "$VCF_LIST")" -eq "$EXPECTED_SAMPLES" ]]

joint_ok=0
if [[ -s "$JOINT" ]] && "$BCFTOOLS" view -h "$JOINT" >/dev/null 2>&1; then
    "$BCFTOOLS" query -l "$JOINT" > "$TMP/joint.samples"
    if cmp -s "$SAMPLE_ORDER" "$TMP/joint.samples" \
        && [[ "$($BCFTOOLS view -H "$JOINT" | wc -l)" -eq "$EXPECTED_RECORDS" ]]; then
        joint_ok=1
    fi
fi

if [[ "$joint_ok" -eq 0 ]]; then
    rm -f "$JOINT" "$JOINT.csi" "$TMP/joint.samples"
    mapfile -t vcfs < "$VCF_LIST"
    "$SMOOVE" paste --name medicago143 --outdir "$COHORT" "${vcfs[@]}"
    generated="$COHORT/medicago143.smoove.square.vcf.gz"
    test -s "$generated"
fi

"$BCFTOOLS" query -l "$JOINT" > "$TMP/joint.samples"
cmp "$SAMPLE_ORDER" "$TMP/joint.samples"
[[ "$($BCFTOOLS view -H "$JOINT" | wc -l)" -eq "$EXPECTED_RECORDS" ]]
"$BCFTOOLS" index -f -c "$JOINT"

if [[ ! -s "$QC" ]] || ! "$BCFTOOLS" view -h "$QC" >/dev/null 2>&1; then
    rm -f "$QC" "$QC.csi"
    "$BCFTOOLS" +setGT "$JOINT" -Ou -- \
        -t q -n . -i 'FMT/GQ<20 || FMT/DP<5' \
    | "$BCFTOOLS" +fill-tags -Ou -- -t AN,AC,AF,MAF,NS \
    | "$BCFTOOLS" view \
        -i 'NS>=115 && MAF>=0.05 && AC>=5 && (AN-AC)>=5 && (INFO/SVTYPE!="BND" || ID~"_1$")' \
        -Oz -o "$QC"
fi

"$BCFTOOLS" index -f -c "$QC"
"$BCFTOOLS" query -l "$QC" > "$TMP/qc.samples"
cmp "$SAMPLE_ORDER" "$TMP/qc.samples"
qc_records="$($BCFTOOLS index -n "$QC")"
[[ "$qc_records" -gt 0 ]]

{
    printf 'svtype\tvariant_count\n'
    "$BCFTOOLS" query -f '%INFO/SVTYPE\n' "$QC" \
        | sort | uniq -c \
        | awk '{print $2 "\t" $1}'
} > "$SUMMARY/sv_cohort_qc_svtype_counts.tsv"

{
    printf 'metric\tvalue\n'
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
    printf 'samples\t%s\n' "$EXPECTED_SAMPLES"
    printf 'input_records\t%s\n' "$EXPECTED_RECORDS"
    printf 'qc_records\t%s\n' "$qc_records"
    printf 'genotype_GQ_min\t20\n'
    printf 'genotype_DP_min\t5\n'
    printf 'site_call_rate_min\t0.80\n'
    printf 'site_called_samples_min\t%s\n' "$MIN_CALLED_SAMPLES"
    printf 'site_MAF_min\t0.05\n'
    printf 'site_MAC_min\t5\n'
    printf 'BND_representation\tone_record_per_MATEID_pair_keep_ID_suffix_1\n'
    printf 'joint_vcf\t%s\n' "$JOINT"
    printf 'qc_vcf\t%s\n' "$QC"
} > "$SUMMARY/sv_cohort_qc_counts.tsv"

cp "$SUMMARY/sv_cohort_qc_counts.tsv" "$SUMMARY/sv_cohort_qc.done"
echo "[$(date '+%F %T')] SV cohort paste and QC completed: records=$qc_records"
