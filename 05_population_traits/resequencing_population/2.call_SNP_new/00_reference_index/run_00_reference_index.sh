#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/38.medicago_resequence
OUT=${BASE}/2.call_SNP_new
REF_SRC=${BASE}/0.ref_genome/ref_genome.T2T.fa
REF_DIR=${OUT}/00_reference_index
REF=${REF_DIR}/ref_genome.T2T.fa
LOG_DIR=${REF_DIR}/logs

CONDA_SH=path/to/home/anaconda3/etc/profile.d/conda.sh
CALL_ENV=call_snp
SAMTOOLS=path/to/software/samtools/samtools-1.15.1/samtools

mkdir -p "${REF_DIR}" "${LOG_DIR}"
cd "${REF_DIR}"

set +u
source "${CONDA_SH}"
conda activate "${CALL_ENV}"
set -u

ln -sfn "${REF_SRC}" "${REF}"

{
  echo "Started: $(date)"
  echo "Host: $(hostname)"
  echo "Reference: ${REF}"
  gatk --version | head -5 || true
  bwa 2>&1 | head -5 || true
  "${SAMTOOLS}" --version | head -3 || true
} > "${LOG_DIR}/00_versions.log" 2>&1

if [[ ! -s "${REF}.fai" ]]; then
  "${SAMTOOLS}" faidx "${REF}" > "${LOG_DIR}/samtools_faidx.log" 2>&1
fi

if [[ ! -s "${REF_DIR}/ref_genome.T2T.dict" ]]; then
  gatk CreateSequenceDictionary \
    -R "${REF}" \
    -O "${REF_DIR}/ref_genome.T2T.dict" \
    > "${LOG_DIR}/gatk_CreateSequenceDictionary.log" 2>&1
fi

if [[ ! -s "${REF}.bwt" ]]; then
  bwa index "${REF}" > "${LOG_DIR}/bwa_index.log" 2>&1
fi

echo "Finished: $(date)" >> "${LOG_DIR}/00_versions.log"
