#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/38.medicago_resequence
OUT=${BASE}/2.call_SNP_new
LOG_DIR=${OUT}/99_run_all_stepwise/logs

mkdir -p "${LOG_DIR}"

{
  echo "Stepwise pipeline started: $(date)"
  echo "Host: $(hostname)"
  echo "OUT=${OUT}"
} > "${LOG_DIR}/run_all_stepwise.progress.log"

bash "${OUT}/02_fastp_all/run.sh" >> "${LOG_DIR}/run_all_stepwise.progress.log" 2>&1
bash "${OUT}/03_bwa_sort_all/run.sh" >> "${LOG_DIR}/run_all_stepwise.progress.log" 2>&1
bash "${OUT}/04_markdup_all/run.sh" >> "${LOG_DIR}/run_all_stepwise.progress.log" 2>&1
bash "${OUT}/05_haplotypecaller_gvcf/run.sh" >> "${LOG_DIR}/run_all_stepwise.progress.log" 2>&1
bash "${OUT}/06_joint_genotype/run.sh" >> "${LOG_DIR}/run_all_stepwise.progress.log" 2>&1
bash "${OUT}/07_filter_snp/run.sh" >> "${LOG_DIR}/run_all_stepwise.progress.log" 2>&1

echo "Stepwise pipeline finished: $(date)" >> "${LOG_DIR}/run_all_stepwise.progress.log"
