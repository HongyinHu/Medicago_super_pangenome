#!/usr/bin/env bash
set -euo pipefail
RUN=path/to/project/N_4.pod_spiny/08_GWAS_pod_spine_SNP_SV_20260706
PHENO=$RUN/input/pod_spine_plink.pheno
PLINK2=path/to/home/anaconda3/envs/call_snp/bin/plink2
PREFIX=$RUN/results/sativa182_Msa_pod_spine_snp.maf005
mkdir -p "$RUN/results" "$RUN/logs" "$RUN/qc"
cd "$RUN"
{
  echo "[START_ASSOC] $(date) host=$(hostname)"
  echo "[STEP] Binary association without covariates"
  "$PLINK2" \
    --pfile "$PREFIX" \
    --allow-extra-chr \
    --pheno "$PHENO" \
    --glm hide-covar \
    --out "$RUN/results/sativa182_pod_spine_snp_gwas.no_covar"
  echo "[STEP] Binary association with PC1-PC5 covariates"
  "$PLINK2" \
    --pfile "$PREFIX" \
    --allow-extra-chr \
    --pheno "$PHENO" \
    --covar "$RUN/qc/sativa182_Msa_pod_spine_snp.maf005.pca.eigenvec" \
    --covar-name PC1,PC2,PC3,PC4,PC5 \
    --glm hide-covar \
    --out "$RUN/results/sativa182_pod_spine_snp_gwas.PC5"
  echo "[STEP] Summarize Chr23997 window Chr4:88,000,000-90,000,000"
  python3 "$RUN/scripts/summarize_plink2_gwas.py" \
    "$RUN/results/sativa182_pod_spine_snp_gwas.no_covar.PHENO1.glm.logistic" \
    Chr4 88000000 90000000 \
    "$RUN/results/Chr23997_window.no_covar.tsv" \
    "$RUN/results/Chr23997_window.no_covar.summary.txt"
  python3 "$RUN/scripts/summarize_plink2_gwas.py" \
    "$RUN/results/sativa182_pod_spine_snp_gwas.PC5.PHENO1.glm.logistic" \
    Chr4 88000000 90000000 \
    "$RUN/results/Chr23997_window.PC5.tsv" \
    "$RUN/results/Chr23997_window.PC5.summary.txt"
  echo "[DONE_ASSOC] $(date)"
  touch "$RUN/results/snp_gwas_plink2.assoc.done"
} 2>&1 | tee "$RUN/logs/run_snp_gwas_plink2_assoc_only.log"
