#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/medicago_resequence_call_SNP_yue_backup
REF=${BASE}/00_reference_index/ref_genome.T2T.fa
RAW_DIR=${BASE}/06_joint_genotype/raw_by_interval
STEP_DIR=${BASE}/08_filter_indel
LOG_DIR=${STEP_DIR}/logs
ENV_BIN=${BASE}/tools/lz_call_snp_env/bin
GATK=${ENV_BIN}/gatk
BCFTOOLS=/usr/local/bin/bcftools

JAVA_MEM=${JAVA_MEM:-24g}
MISS_FRACTION=${MISS_FRACTION:-0.2}
MIN_DP=${MIN_DP:-6}
MAX_DP=${MAX_DP:-85}
PARALLEL_FILTER=${PARALLEL_FILTER:-8}
RUN_ID=${RUN_ID:-$(date +%Y%m%d_%H%M%S)}
WORK_DIR=${STEP_DIR}/by_interval_${RUN_ID}
FINAL=${STEP_DIR}/DP_${MIN_DP}-${MAX_DP}_miss_${MISS_FRACTION}.all_samples.biallelic.INDEL.vcf.gz
ORDER=${WORK_DIR}/final_order.list

export PATH=${ENV_BIN}:${PATH}
export JAVA_HOME
JAVA_HOME=$(dirname "$(dirname "$(readlink -f "${ENV_BIN}/java")")")

mkdir -p "${LOG_DIR}" "${WORK_DIR}/parts" "${WORK_DIR}/logs"
cd "${BASE}"

for f in "${REF}" "${RAW_DIR}/Chr1.raw.vcf.gz" "${RAW_DIR}/Chr8.raw.vcf.gz"; do
  [[ -s "$f" ]] || { echo "Missing required file: $f" >&2; exit 1; }
done
for f in "${GATK}" "${BCFTOOLS}"; do
  [[ -x "$f" ]] || { echo "Missing executable: $f" >&2; exit 1; }
done

contigs=(Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 scaffold1 scaffold2 scaffold3 scaffold4 scaffold5 scaffold6 scaffold7)
: > "${ORDER}"
for c in "${contigs[@]}"; do
  raw="${RAW_DIR}/${c}.raw.vcf.gz"
  [[ -s "$raw" && -s "$raw.tbi" ]] || { echo "Missing raw/tbi for $c: $raw" >&2; exit 1; }
  echo "${WORK_DIR}/parts/${c}.filtered.biallelic.INDEL.vcf.gz" >> "${ORDER}"
done

log_step() { echo "[$(date '+%F %T %Z')] $*"; }
check_vcf() {
  local f="$1"
  [[ -s "$f" && -s "$f.tbi" ]] || { echo "Missing output or tbi: $f" >&2; exit 1; }
}
index_vcf_if_needed() {
  local f="$1"
  [[ -s "$f.tbi" ]] || "${BCFTOOLS}" index -t "$f"
}

filter_one() {
  local c="$1"
  local raw="${RAW_DIR}/${c}.raw.vcf.gz"
  local prefix="${WORK_DIR}/parts/${c}"
  local log="${WORK_DIR}/logs/${c}.log"
  local indel_raw="${prefix}.raw.INDEL.vcf.gz"
  local hard="${prefix}.raw.INDEL.hard_filter.vcf.gz"
  local pass="${prefix}.raw.INDEL.pass.vcf.gz"
  local gtdp="${prefix}.raw.INDEL.pass.GT_DP_${MIN_DP}_${MAX_DP}.vcf.gz"
  local final="${prefix}.filtered.biallelic.INDEL.vcf.gz"

  {
    log_step "START ${c}"
    if [[ ! -s "$final" || ! -s "$final.tbi" ]]; then
      "${GATK}" --java-options "-Xmx${JAVA_MEM}" SelectVariants \
        -R "${REF}" -V "$raw" \
        --select-type-to-include INDEL \
        --restrict-alleles-to BIALLELIC \
        -O "$indel_raw"
      index_vcf_if_needed "$indel_raw"
      check_vcf "$indel_raw"

      "${GATK}" --java-options "-Xmx${JAVA_MEM}" VariantFiltration \
        -R "${REF}" -V "$indel_raw" -O "$hard" \
        --filter-name "QD_lt_2" --filter-expression "QD < 2.0" \
        --filter-name "QUAL_lt_30" --filter-expression "QUAL < 30.0" \
        --filter-name "FS_gt_200" --filter-expression "FS > 200.0" \
        --filter-name "ReadPosRankSum_lt_-20" --filter-expression "ReadPosRankSum < -20.0"
      index_vcf_if_needed "$hard"
      check_vcf "$hard"
      rm -f "$indel_raw" "$indel_raw.tbi"

      "${GATK}" --java-options "-Xmx${JAVA_MEM}" SelectVariants \
        -R "${REF}" -V "$hard" --exclude-filtered true -O "$pass"
      index_vcf_if_needed "$pass"
      check_vcf "$pass"
      rm -f "$hard" "$hard.tbi"

      "${GATK}" --java-options "-Xmx${JAVA_MEM}" VariantFiltration \
        -R "${REF}" -V "$pass" -O "$gtdp" \
        --genotype-filter-name "DP_outside_${MIN_DP}_${MAX_DP}" \
        --genotype-filter-expression "DP < ${MIN_DP} || DP > ${MAX_DP}" \
        --set-filtered-genotype-to-no-call
      index_vcf_if_needed "$gtdp"
      check_vcf "$gtdp"
      rm -f "$pass" "$pass.tbi"

      "${GATK}" --java-options "-Xmx${JAVA_MEM}" SelectVariants \
        -R "${REF}" -V "$gtdp" --max-nocall-fraction "${MISS_FRACTION}" -O "$final"
      index_vcf_if_needed "$final"
      check_vcf "$final"
      rm -f "$gtdp" "$gtdp.tbi"
    fi
    check_vcf "$final"
    log_step "FINISHED ${c}"
  } > "$log" 2>&1
}

export -f filter_one log_step check_vcf index_vcf_if_needed
export BASE REF RAW_DIR STEP_DIR LOG_DIR GATK BCFTOOLS JAVA_MEM MISS_FRACTION MIN_DP MAX_DP WORK_DIR

log_step "START yue by-interval INDEL filter RUN_ID=${RUN_ID} PARALLEL_FILTER=${PARALLEL_FILTER} JAVA_MEM=${JAVA_MEM}"
printf '%s\n' "${contigs[@]}" | xargs -n 1 -P "${PARALLEL_FILTER}" bash -c 'filter_one "$0"'

while read -r f; do check_vcf "$f"; done < "${ORDER}"
log_step "Concat filtered interval INDEL VCFs"
"${BCFTOOLS}" concat -f "${ORDER}" -Oz -o "${FINAL}"
"${BCFTOOLS}" index -t "${FINAL}"
check_vcf "${FINAL}"
"${BCFTOOLS}" stats "${FINAL}" > "${STEP_DIR}/DP_${MIN_DP}-${MAX_DP}_miss_${MISS_FRACTION}.all_samples.biallelic.INDEL.stats.txt"
ln -sfn "${FINAL}" "${BASE}/DP_${MIN_DP}-${MAX_DP}_miss_${MISS_FRACTION}.all_samples.biallelic.INDEL.vcf.gz"
ln -sfn "${FINAL}.tbi" "${BASE}/DP_${MIN_DP}-${MAX_DP}_miss_${MISS_FRACTION}.all_samples.biallelic.INDEL.vcf.gz.tbi"

# Per-contig INDEL files are temporary shards. Remove only after final VCF and index exist.
if [[ -s "${FINAL}" && -s "${FINAL}.tbi" ]]; then
  rm -f "${WORK_DIR}/parts/"*.filtered.biallelic.INDEL.vcf.gz "${WORK_DIR}/parts/"*.filtered.biallelic.INDEL.vcf.gz.tbi
fi
log_step "FINISHED yue by-interval INDEL filter: ${FINAL}"
