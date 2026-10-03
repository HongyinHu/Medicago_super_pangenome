#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/03_population_structure"
QC="$ROOT/02_snp_indel_qc/results/plink"
PHENO="$ROOT/01_input_audit/summary/phenotype_144.tsv"
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma
PYTHON=/usr/bin/python3

mkdir -p "$STAGE/results/plink" "$STAGE/results/gemma" "$STAGE/summary" "$STAGE/logs" "$STAGE/tmp"

SNP_IN="$QC/144.SNP.common"
INDEL_IN="$QC/144.INDEL.common"
SNP="$STAGE/results/plink/analysis.SNP"
INDEL="$STAGE/results/plink/analysis.INDEL"
LD="$STAGE/results/plink/analysis.SNP.ldpruned"
PRUNE="$STAGE/results/plink/SNP.ld"
PCA="$STAGE/results/plink/SNP.pca"

for suffix in bed bim fam; do
  test -s "$SNP_IN.$suffix"
  test -s "$INDEL_IN.$suffix"
done

awk '{print $1, $2}' "$SNP_IN.fam" > "$STAGE/tmp/snp.samples"
awk '{print $1, $2}' "$INDEL_IN.fam" > "$STAGE/tmp/indel.samples"
if ! cmp -s "$STAGE/tmp/snp.samples" "$STAGE/tmp/indel.samples"; then
  echo "SNP and INDEL sample order differs" >&2
  diff -u "$STAGE/tmp/snp.samples" "$STAGE/tmp/indel.samples" >&2 || true
  exit 2
fi

if [[ ! -s "$SNP.bed" ]]; then
  "$PLINK" --bfile "$SNP_IN" --chr-set 8 no-xy --make-bed --out "$SNP"
fi
if [[ ! -s "$INDEL.bed" ]]; then
  "$PLINK" --bfile "$INDEL_IN" --chr-set 8 no-xy --make-bed --out "$INDEL"
fi

if [[ ! -s "$PRUNE.prune.in" ]]; then
  echo "[$(date '+%F %T')] LD pruning SNPs"
  "$PLINK" --bfile "$SNP" --chr-set 8 no-xy --indep-pairwise 50 5 0.2 --out "$PRUNE"
fi
if [[ ! -s "$LD.bed" ]]; then
  "$PLINK" --bfile "$SNP" --chr-set 8 no-xy --extract "$PRUNE.prune.in" --make-bed --out "$LD"
fi
if [[ ! -s "$PCA.eigenvec" ]]; then
  echo "[$(date '+%F %T')] Calculating SNP PCA"
  "$PLINK" --bfile "$LD" --chr-set 8 no-xy --pca 10 --out "$PCA"
fi

"$PYTHON" "$STAGE/scripts/build_structure_inputs.py" \
  "$SNP.fam" "$PHENO" "$PCA.eigenvec" "$STAGE/summary" 5

KIN="$STAGE/results/gemma/output/SNP_kinship.cXX.txt"
if [[ ! -s "$KIN" ]]; then
  echo "[$(date '+%F %T')] Calculating centered SNP kinship with GEMMA"
  (
    cd "$STAGE/results/gemma"
    "$GEMMA" -bfile "$LD" -gk 1 -o SNP_kinship
  )
fi
test -s "$KIN"

{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'samples\t%s\n' "$(wc -l < "$SNP.fam")"
  printf 'common_snp_variants\t%s\n' "$(wc -l < "$SNP.bim")"
  printf 'common_indel_variants\t%s\n' "$(wc -l < "$INDEL.bim")"
  printf 'ld_pruned_snp_variants\t%s\n' "$(wc -l < "$LD.bim")"
  printf 'kinship\t%s\n' "$KIN"
  printf 'covariates\t%s\n' "$STAGE/summary/covariates.pc5.txt"
} > "$STAGE/summary/population_structure.done"

echo "[$(date '+%F %T')] Population structure completed"
