#!/usr/bin/env bash
set -euo pipefail

LABEL=${1:?Usage: run_marker_gwas.sh SNP_or_INDEL}
[[ "$LABEL" == SNP || "$LABEL" == INDEL ]] || { echo "LABEL must be SNP or INDEL" >&2; exit 2; }
CHR=${SLURM_ARRAY_TASK_ID:-${2:-}}
[[ "$CHR" =~ ^[1-8]$ ]] || { echo "Chromosome must be supplied as array task 1-8 or second argument" >&2; exit 2; }

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript
export GWAS_ROOT="$ROOT"
export R_LIBS_USER="$ROOT/software/Rlib"
export TMPDIR="$ROOT/tmp/R_${LABEL}_chr${CHR}"
mkdir -p "$TMPDIR" "$ROOT/02_snp_indel_gwas/gmmat" "$ROOT/02_snp_indel_gwas/logistic" "$ROOT/logs" "$ROOT/summary"
LOG="$ROOT/logs/run_${LABEL}_chr${CHR}_gmmat.log"
exec > >(tee -a "$LOG") 2>&1

[[ -s "$ROOT/summary/structure_144.done" ]]
BFILE="$ROOT/01_structure/plink/coil144.${LABEL}.analysis.chr${CHR}"
MODEL="$ROOT/01_structure/model/${LABEL}.model.tsv"
GRM="$ROOT/01_structure/plink/coil144.SNP.grm"
OUT="$ROOT/02_snp_indel_gwas/gmmat/${LABEL}.chr${CHR}.coil_ge1.pc5.gmmat.score.tsv"

echo "[$(date '+%F %T %Z')] $LABEL Chr${CHR} GMMAT start on $(hostname)"
if [[ ! -s "$OUT" ]]; then
  "$RSCRIPT" --vanilla "$ROOT/scripts/run_gmmat_gwas.R" \
    "$BFILE" "$MODEL" "$GRM.rel" "$GRM.rel.id" "$OUT" "$LABEL" 5 "${SLURM_CPUS_PER_TASK:-8}"
fi
test -s "$OUT"

SENS="$ROOT/02_snp_indel_gwas/logistic/${LABEL}.chr${CHR}.coil_ge1.pc5"
if [[ ! -s "$SENS.assoc.logistic" ]]; then
  "$PLINK" --bfile "$BFILE" --chr-set 8 no-xy --allow-no-sex \
    --pheno "$ROOT/01_structure/model/${LABEL}.plink.pheno.tsv" \
    --covar "$ROOT/01_structure/model/${LABEL}.plink.covar.tsv" \
    --covar-name PC1,PC2,PC3,PC4,PC5 --logistic hide-covar --ci 0.95 --out "$SENS"
fi

{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'variant_class\t%s\n' "$LABEL"
  printf 'primary_model\tGMMAT binomial logistic mixed model; SNP GRM; PC1-PC5\n'
  printf 'sensitivity_model\tPLINK ordinary logistic regression; PC1-PC5\n'
  printf 'samples\t%s\n' "$(wc -l < "$BFILE.fam")"
  printf 'markers\t%s\n' "$(wc -l < "$BFILE.bim")"
  printf 'chromosome\t%s\n' "$CHR"
} > "$ROOT/summary/${LABEL}.chr${CHR}.done"
echo "[$(date '+%F %T %Z')] $LABEL Chr${CHR} GMMAT complete"
