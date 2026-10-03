#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
SENS=$ROOT/05_model_sensitivity_20260715
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript
export GWAS_ROOT="$ROOT"
export R_LIBS_USER="$ROOT/software/Rlib"
export TMPDIR="$SENS/tmp/SV"

mkdir -p "$TMPDIR" "$SENS/model_B_plink_pc3/SV" \
  "$SENS/model_C_gmmat_grm_pc2/SV" "$SENS/logs" "$SENS/summary/tasks"
LOG="$SENS/logs/SV.model_sensitivity.log"
exec > >(tee -a "$LOG") 2>&1

test -s "$ROOT/summary/SV_gwas.done"
BFILE="$ROOT/03_sv_gwas/plink/coil144.SV.GQ20.analysis"
MODEL="$ROOT/03_sv_gwas/model/SV.model.tsv"
PHENO="$ROOT/03_sv_gwas/model/SV.plink.pheno.tsv"
COVAR="$ROOT/03_sv_gwas/model/SV.plink.covar.tsv"
GRM="$ROOT/01_structure/plink/coil144.SNP.grm"
test -s "$BFILE.bed"

echo "[$(date '+%F %T %Z')] SV Model B start on $(hostname)"
BOUT="$SENS/model_B_plink_pc3/SV/SV.pc3"
if [[ ! -s "$BOUT.assoc.logistic" ]]; then
  "$PLINK" --bfile "$BFILE" --chr-set 8 no-xy --allow-no-sex \
    --pheno "$PHENO" --covar "$COVAR" --covar-name PC1,PC2,PC3 \
    --logistic hide-covar --ci 0.95 --out "$BOUT"
fi
test -s "$BOUT.assoc.logistic"

echo "[$(date '+%F %T %Z')] SV Model C start"
COUT="$SENS/model_C_gmmat_grm_pc2/SV/SV.grm_pc2.gmmat.score.tsv"
if [[ ! -s "$COUT" ]]; then
  "$RSCRIPT" --vanilla "$ROOT/scripts/run_gmmat_gwas.R" \
    "$BFILE" "$MODEL" "$GRM.rel" "$GRM.rel.id" "$COUT" \
    SV_Model_C 2 "${SLURM_CPUS_PER_TASK:-2}"
fi
test -s "$COUT"

{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'model_B\tPLINK logistic; PC1-PC3; no GRM\n'
  printf 'model_C\tGMMAT binomial mixed model; PC1-PC2; SNP GRM\n'
  printf 'samples\t%s\n' "$(wc -l < "$BFILE.fam")"
  printf 'markers\t%s\n' "$(wc -l < "$BFILE.bim")"
} > "$SENS/summary/tasks/SV.done"
echo "[$(date '+%F %T %Z')] SV sensitivity complete"
