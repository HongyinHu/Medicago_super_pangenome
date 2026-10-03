#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/07_joint_gwas"
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
RSCRIPT=path/to/home/anaconda3/envs/panpop/bin/Rscript
FAI=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa.fai
SNP_INDEL="$ROOT/04_snp_indel_gwas/results/summary/SNP_INDEL.primary.assoc.tsv.gz"
SV="$ROOT/06_sv_gwas/results/summary/SV.primary.assoc.tsv.gz"

mkdir -p "$STAGE/results/tables" "$STAGE/results/figures" "$STAGE/logs/slurm" "$STAGE/summary"
test -s "$ROOT/04_snp_indel_gwas/summary/snp_indel_gwas.done"
test -s "$ROOT/06_sv_gwas/summary/sv_gwas.done"
test -s "$SNP_INDEL"
test -s "$SV"
test -s "$FAI"

"$PYTHON" "$STAGE/scripts/integrate_gwas_results.py" \
    --snp-indel "$SNP_INDEL" \
    --sv "$SV" \
    --outdir "$STAGE/results/tables"

"$RSCRIPT" "$STAGE/scripts/plot_joint_gwas.R" \
    "$STAGE/results/tables/SNP_INDEL_SV.primary.assoc.tsv.gz" \
    "$FAI" \
    "$STAGE/results/figures"

test -s "$STAGE/results/tables/joint_association_summary.tsv"
test -s "$STAGE/results/figures/SNP_INDEL_SV.manhattan.png"
test -s "$STAGE/results/figures/SNP_INDEL_SV.QQ.png"
test -s "$STAGE/results/figures/Chr23997.plusminus1Mb.png"

{
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
    printf 'model\tGEMMA_LMM_SNP_kinship_PC1-PC5\n'
    printf 'samples\t143\n'
    printf 'joint_table\t%s\n' "$STAGE/results/tables/SNP_INDEL_SV.primary.assoc.tsv.gz"
    printf 'figures\t%s\n' "$STAGE/results/figures"
} > "$STAGE/summary/integration.done"

echo "[$(date '+%F %T')] Joint GWAS integration completed"
