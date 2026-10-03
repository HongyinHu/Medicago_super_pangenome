#!/usr/bin/env bash
set -euo pipefail

RUN="$1"
TASK_LINE="$2"
THREADS="$3"

MINIMAP2="path/to/home/anaconda3/envs/minimap2_env/bin/minimap2"
SAMTOOLS="path/to/home/anaconda3/envs/biosofeware/bin/samtools"
SNIFFLES2="path/to/home/anaconda3/envs/sniffles_env/bin/sniffles"
CUTESV="path/to/home/anaconda3/envs/cuteSV/bin/cuteSV"
PBSV="path/to/home/anaconda3/envs/biosofeware/bin/pbsv"
BGZIP="path/to/home/anaconda3/envs/biosofeware/bin/bgzip"
TABIX="path/to/home/anaconda3/envs/biosofeware/bin/tabix"

IFS=$'\t' read -r TASK_ID SAMPLE ORIGINAL REF_ID REF_FA REF_MMI READS ROLE COILING SPINE READ_GB <<< "${TASK_LINE}"

HOST=$(hostname)
TASK_LOG="${RUN}/logs/tasks/${TASK_ID}.${HOST}.$(date +%Y%m%d_%H%M%S).log"
STATUS_DIR="${RUN}/status"
TMP_ROOT="${RUN}/tmp/${HOST}/${TASK_ID}"

mkdir -p "${TMP_ROOT}" "${STATUS_DIR}" \
  "${RUN}/03_per_sample/${REF_ID}/bam" \
  "${RUN}/03_per_sample/${REF_ID}/sniffles2" \
  "${RUN}/03_per_sample/${REF_ID}/cutesv" \
  "${RUN}/03_per_sample/${REF_ID}/pbsv" \
  "${RUN}/03_per_sample/${REF_ID}/snf"

exec > >(tee -a "${TASK_LOG}") 2>&1

fail_task() {
  local code=$?
  echo "FAILED task=${TASK_ID} code=${code} host=${HOST} time=$(date -Is)"
  echo -e "${TASK_ID}\t${HOST}\t${code}\t$(date -Is)\t${TASK_LOG}" > "${STATUS_DIR}/${TASK_ID}.failed"
  exit "${code}"
}
trap fail_task ERR

echo "START task=${TASK_ID} sample=${SAMPLE} original=${ORIGINAL} ref=${REF_ID} host=${HOST} threads=${THREADS} time=$(date -Is)"
echo "reads=${READS}"

if [[ ! -s "${READS}" ]]; then
  echo "Missing reads: ${READS}" >&2
  echo -e "${TASK_ID}\t${HOST}\t2\t$(date -Is)\t${TASK_LOG}" > "${STATUS_DIR}/${TASK_ID}.failed"
  exit 2
fi
if [[ ! -s "${REF_FA}" || ! -s "${REF_MMI}" ]]; then
  echo "Missing reference/index: ${REF_FA} ${REF_MMI}" >&2
  echo -e "${TASK_ID}\t${HOST}\t3\t$(date -Is)\t${TASK_LOG}" > "${STATUS_DIR}/${TASK_ID}.failed"
  exit 3
fi

THREADS=$(( THREADS < 2 ? 2 : THREADS ))
SORT_THREADS=4
if (( THREADS < 10 )); then
  SORT_THREADS=2
fi
MAP_THREADS=$(( THREADS - SORT_THREADS ))
if (( MAP_THREADS < 1 )); then
  MAP_THREADS=1
fi
CALL_THREADS=$(( THREADS / 2 ))
if (( CALL_THREADS < 2 )); then
  CALL_THREADS=2
fi
if (( CALL_THREADS > 12 )); then
  CALL_THREADS=12
fi

BAM="${RUN}/03_per_sample/${REF_ID}/bam/${SAMPLE}.sorted.bam"
SNIFFLES_VCF="${RUN}/03_per_sample/${REF_ID}/sniffles2/${SAMPLE}.sniffles2.vcf.gz"
SNIFFLES_SNF="${RUN}/03_per_sample/${REF_ID}/snf/${SAMPLE}.sniffles2.snf"
CUTESV_VCF_RAW="${RUN}/03_per_sample/${REF_ID}/cutesv/${SAMPLE}.cutesv.vcf"
CUTESV_VCF="${CUTESV_VCF_RAW}.gz"
CUTESV_WORK="${RUN}/tmp/${HOST}/${TASK_ID}.cutesv_work.$(date +%s).$$"
PBSV_SIG="${RUN}/03_per_sample/${REF_ID}/pbsv/${SAMPLE}.svsig.gz"
PBSV_VCF_RAW="${RUN}/03_per_sample/${REF_ID}/pbsv/${SAMPLE}.pbsv.vcf"
PBSV_VCF="${PBSV_VCF_RAW}.gz"

if [[ ! -s "${BAM}.bai" ]]; then
  echo "MAPPING ${TASK_ID} map_threads=${MAP_THREADS} sort_threads=${SORT_THREADS}"
  "${MINIMAP2}" -ax map-hifi -Y -t "${MAP_THREADS}" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" "${REF_MMI}" "${READS}" \
    | "${SAMTOOLS}" sort -@ "${SORT_THREADS}" -m 2G -T "${TMP_ROOT}/${SAMPLE}" -o "${BAM}" -
  "${SAMTOOLS}" index -@ "${SORT_THREADS}" "${BAM}"
else
  echo "SKIP_MAPPING existing ${BAM}.bai"
fi

if [[ ! -s "${SNIFFLES_VCF}" ]]; then
  echo "SNIFFLES2 ${TASK_ID}"
  "${SNIFFLES2}" --input "${BAM}" --vcf "${SNIFFLES_VCF}" --snf "${SNIFFLES_SNF}" \
    --reference "${REF_FA}" --threads "${CALL_THREADS}" --minsvlen 50 --mapq 20 --allow-overwrite \
    --tmp-dir "${TMP_ROOT}"
  "${TABIX}" -f -p vcf "${SNIFFLES_VCF}" || true
else
  echo "SKIP_SNIFFLES2 existing ${SNIFFLES_VCF}"
fi

if [[ ! -s "${CUTESV_VCF}" ]]; then
  echo "CUTESV ${TASK_ID}"
  mkdir -p "${CUTESV_WORK}"
  "${CUTESV}" "${BAM}" "${REF_FA}" "${CUTESV_VCF_RAW}" "${CUTESV_WORK}" \
    -t "${CALL_THREADS}" -S "${SAMPLE}" --genotype -q 20 -s 3 -l 50 \
    --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
  "${BGZIP}" -f "${CUTESV_VCF_RAW}"
  "${TABIX}" -f -p vcf "${CUTESV_VCF}"
else
  echo "SKIP_CUTESV existing ${CUTESV_VCF}"
fi

if [[ ! -s "${PBSV_VCF}" ]]; then
  echo "PBSV ${TASK_ID}"
  if [[ ! -s "${PBSV_SIG}" ]]; then
    "${PBSV}" discover --hifi -s "${SAMPLE}" "${BAM}" "${PBSV_SIG}"
  fi
  "${PBSV}" call --hifi -j "${CALL_THREADS}" "${REF_FA}" "${PBSV_SIG}" "${PBSV_VCF_RAW}"
  "${BGZIP}" -f "${PBSV_VCF_RAW}"
  "${TABIX}" -f -p vcf "${PBSV_VCF}"
else
  echo "SKIP_PBSV existing ${PBSV_VCF}"
fi

echo -e "${TASK_ID}\t${SAMPLE}\t${ORIGINAL}\t${REF_ID}\t${HOST}\t${THREADS}\t$(date -Is)\t${TASK_LOG}" > "${STATUS_DIR}/${TASK_ID}.done"
rm -f "${STATUS_DIR}/${TASK_ID}.failed"
echo "DONE task=${TASK_ID} host=${HOST} time=$(date -Is)"
