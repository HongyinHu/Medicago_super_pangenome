#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/07_joint_gwas"
test -s "$STAGE/summary/integration.done"
test -s "$STAGE/summary/sensitivity.done"
test -s "$STAGE/results/tables/joint_association_summary.tsv"
test -s "$STAGE/results/tables/model_sensitivity.tsv"
test -s "$STAGE/results/figures/SNP_INDEL_SV.manhattan.png"
test -s "$STAGE/results/figures/SNP_INDEL_SV.QQ.png"
test -s "$STAGE/results/figures/Chr23997.plusminus1Mb.png"

{
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
    printf 'samples\t143\n'
    printf 'main_model\tSNP_kinship_plus_PC1-PC5\n'
    printf 'sensitivity_models\tSNP_kinship_only,SNP_kinship_plus_PC1-PC3\n'
    printf 'Chr23997_exact_221bp_DEL_jointly_genotyped\tno\n'
    printf 'results\t%s\n' "$STAGE/results"
} > "$STAGE/summary/joint_gwas.complete"

echo "[$(date '+%F %T')] Stage 07 complete"
