#!/usr/bin/env bash
set -euo pipefail
RUN=path/to/project/N_4.pod_spiny/08_GWAS_pod_spine_SNP_SV_20260706
PHENO=$RUN/input/pod_spine_plink2_header.pheno
PLINK2=path/to/home/anaconda3/envs/call_snp/bin/plink2
PREFIX=$RUN/results/sativa182_Msa_pod_spine_snp.maf005
cd "$RUN"
{
  echo "[START_CORRECTED] $(date) host=$(hostname)"
  "$PLINK2" --pfile "$PREFIX" --allow-extra-chr --pheno "$PHENO" --pheno-name pod_spine --glm hide-covar --out "$RUN/results/sativa182_pod_spine_snp_gwas.corrected.no_covar"
  "$PLINK2" --pfile "$PREFIX" --allow-extra-chr --pheno "$PHENO" --pheno-name pod_spine --covar "$RUN/qc/sativa182_Msa_pod_spine_snp.maf005.pca.eigenvec" --covar-name PC1,PC2,PC3,PC4,PC5 --glm hide-covar --out "$RUN/results/sativa182_pod_spine_snp_gwas.corrected.PC5"
  python3 "$RUN/scripts/summarize_plink2_gwas.py" "$RUN/results/sativa182_pod_spine_snp_gwas.corrected.no_covar.pod_spine.glm.logistic" Chr4 88000000 90000000 "$RUN/results/Chr23997_window.corrected.no_covar.tsv" "$RUN/results/Chr23997_window.corrected.no_covar.summary.txt"
  python3 "$RUN/scripts/summarize_plink2_gwas.py" "$RUN/results/sativa182_pod_spine_snp_gwas.corrected.PC5.pod_spine.glm.logistic" Chr4 88000000 90000000 "$RUN/results/Chr23997_window.corrected.PC5.tsv" "$RUN/results/Chr23997_window.corrected.PC5.summary.txt"
  echo "[DONE_CORRECTED] $(date)"
  touch "$RUN/results/snp_gwas_plink2.corrected.done"
} 2>&1 | tee "$RUN/logs/run_snp_gwas_plink2_assoc_corrected.log"
