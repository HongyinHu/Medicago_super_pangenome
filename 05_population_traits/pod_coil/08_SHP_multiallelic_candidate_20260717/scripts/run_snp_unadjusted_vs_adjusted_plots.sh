#!/usr/bin/env bash
set -eo pipefail

source ~/.bashrc
conda activate R

RUN=path/to/project/N_5.pod_coil/08_SHP_multiallelic_candidate_20260717
GWAS=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714

Rscript "$RUN/scripts/plot_snp_unadjusted_vs_adjusted.R" \
  "$RUN/uncorrected/SNP_unadjusted.assoc" \
  "$GWAS/02_snp_indel_gwas/gmmat/SNP.coil_ge1.pc5.gmmat.score.tsv" \
  "$RUN/figures/SNP_unadjusted_vs_adjusted_20260717"
