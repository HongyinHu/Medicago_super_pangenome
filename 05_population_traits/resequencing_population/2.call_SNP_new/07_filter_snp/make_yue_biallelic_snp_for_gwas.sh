#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/medicago_resequence_call_SNP_yue_backup
STEP_DIR=${BASE}/07_filter_snp
BCFTOOLS=/usr/local/bin/bcftools
THREADS=${THREADS:-16}
RUN_ID=${RUN_ID:-$(date +%Y%m%d_%H%M%S)}

INPUT=${STEP_DIR}/DP_6-85_miss_0.2.all_samples.SNP.vcf.gz
FINAL=${STEP_DIR}/DP_6-85_miss_0.2.all_samples.biallelic.SNP.vcf.gz
TMP=${STEP_DIR}/.DP_6-85_miss_0.2.all_samples.biallelic.SNP.${RUN_ID}.tmp.vcf.gz
STATS=${STEP_DIR}/DP_6-85_miss_0.2.all_samples.biallelic.SNP.stats.txt

[[ -s "${INPUT}" && -s "${INPUT}.tbi" ]] || { echo "Missing filtered SNP input: ${INPUT}" >&2; exit 1; }
[[ -x "${BCFTOOLS}" ]] || { echo "Missing bcftools: ${BCFTOOLS}" >&2; exit 1; }

rm -f "${TMP}" "${TMP}.tbi"
echo "[$(date '+%F %T %Z')] START biallelic SNP extraction RUN_ID=${RUN_ID} THREADS=${THREADS}"
"${BCFTOOLS}" view --threads "${THREADS}" -m2 -M2 -v snps -Oz -o "${TMP}" "${INPUT}"
"${BCFTOOLS}" index --threads "${THREADS}" -t "${TMP}"

samples=$("${BCFTOOLS}" query -l "${TMP}" | wc -l)
[[ "${samples}" -eq 184 ]] || { echo "Unexpected sample count: ${samples}" >&2; exit 2; }
records=$("${BCFTOOLS}" index -n "${TMP}")
[[ "${records}" -gt 0 ]] || { echo "No biallelic SNP records produced" >&2; exit 3; }

mv -f "${TMP}" "${FINAL}"
mv -f "${TMP}.tbi" "${FINAL}.tbi"
"${BCFTOOLS}" stats --threads "${THREADS}" "${FINAL}" > "${STATS}"
ln -sfn "${FINAL}" "${BASE}/DP_6-85_miss_0.2.all_samples.biallelic.SNP.vcf.gz"
ln -sfn "${FINAL}.tbi" "${BASE}/DP_6-85_miss_0.2.all_samples.biallelic.SNP.vcf.gz.tbi"
echo "[$(date '+%F %T %Z')] FINISHED biallelic SNP extraction samples=${samples} records=${records} output=${FINAL}"
