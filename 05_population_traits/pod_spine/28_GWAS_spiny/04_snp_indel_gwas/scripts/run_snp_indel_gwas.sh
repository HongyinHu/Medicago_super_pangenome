#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/04_snp_indel_gwas"
STRUCT="$ROOT/03_population_structure"
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python

mkdir -p "$STAGE/results/gemma" "$STAGE/results/summary" "$STAGE/logs" "$STAGE/summary"

test -s "$STRUCT/summary/population_structure.done"
KIN="$STRUCT/results/gemma/output/SNP_kinship.cXX.txt"
COVAR="$STRUCT/summary/covariates.pc5.txt"
PHENO="$STRUCT/summary/phenotype.0_1.txt"
test -s "$KIN"
test -s "$COVAR"
test -s "$PHENO"

run_one() {
  local label=$1
  local prefix="$STRUCT/results/plink/analysis.${label}"
  local result="$STAGE/results/gemma/output/${label}.pc5.assoc.txt"
  if [[ ! -s "$result" ]]; then
    echo "[$(date '+%F %T')] GEMMA ${label} LMM start"
    (
      cd "$STAGE/results/gemma"
      "$GEMMA" -bfile "$prefix" -k "$KIN" -c "$COVAR" -p "$PHENO" -lmm 4 -o "${label}.pc5"
    )
  fi
  test -s "$result"
}

run_one SNP
run_one INDEL

"$PYTHON" "$STAGE/scripts/summarize_gemma.py" \
  "$STAGE/results/gemma/output/SNP.pc5.assoc.txt" \
  "$STAGE/results/gemma/output/INDEL.pc5.assoc.txt" \
  "$STAGE/results/summary"

{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'samples\t%s\n' "$(wc -l < "$STRUCT/results/plink/analysis.SNP.fam")"
  printf 'phenotype\t0=spineless,1=spiny\n'
  printf 'model\tGEMMA_LMM_Wald_primary\n'
  printf 'kinship\tSNP_LD_pruned_centered\n'
  printf 'covariates\tintercept_plus_PC1-PC5\n'
} > "$STAGE/summary/snp_indel_gwas.done"

echo "[$(date '+%F %T')] SNP/INDEL GWAS completed"
