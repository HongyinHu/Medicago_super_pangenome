#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/02_SNP_INDEL_SV_GWAS_20260713
SOURCE=path/to/project/N_4.pod_spiny/28_GWAS_spiny
PHENO=path/to/project/N_4.pod_spiny/21_resequence_phenotype_completed_20260708/results/resequence_sample_pod_traits_completed.tsv
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python

mkdir -p "$ROOT"/{00_input,01_structure,02_snp_indel_gwas,03_sv_gwas,04_candidate_gene_peaks,logs,summary,tmp}
mkdir -p "$ROOT/01_structure"/{plink,gemma,summary} "$ROOT/02_snp_indel_gwas"/{gemma,results}
LOG="$ROOT/logs/run_snp_indel_gwas.log"
exec > >(tee -a "$LOG") 2>&1

echo "[$(date '+%F %T %Z')] strict coil SNP/INDEL GWAS start on $(hostname)"

SOURCE_SNP="$SOURCE/03_population_structure/results/plink/analysis.SNP"
SOURCE_INDEL="$SOURCE/03_population_structure/results/plink/analysis.INDEL"
SOURCE_LD="$SOURCE/03_population_structure/results/plink/analysis.SNP.ldpruned"

"$PYTHON" "$ROOT/scripts/prepare_strict_coil.py" \
  "$PHENO" "$SOURCE_SNP.fam" "$ROOT/00_input"

KEEP="$ROOT/00_input/strict_samples.keep"
STRICT_PHENO="$ROOT/00_input/coil_phenotype_strict.tsv"

for label in SNP INDEL; do
  src="$SOURCE/03_population_structure/results/plink/analysis.$label"
  out="$ROOT/01_structure/plink/coil.$label"
  if [[ ! -s "$out.bed" ]]; then
    "$PLINK" --bfile "$src" --keep "$KEEP" --chr-set 8 no-xy \
      --mind 0.20 --geno 0.20 --maf 0.05 --make-bed --out "$out"
  fi
done

LD="$ROOT/01_structure/plink/coil.SNP.ldpruned"
if [[ ! -s "$LD.bed" ]]; then
  "$PLINK" --bfile "$SOURCE_LD" --keep "$KEEP" --chr-set 8 no-xy \
    --mind 0.20 --geno 0.20 --maf 0.05 --make-bed --out "$LD"
fi

PCA="$ROOT/01_structure/plink/coil.SNP.pca"
if [[ ! -s "$PCA.eigenvec" ]]; then
  "$PLINK" --bfile "$LD" --chr-set 8 no-xy --pca 10 --out "$PCA"
fi

"$PYTHON" "$ROOT/scripts/build_gemma_inputs.py" \
  "$ROOT/01_structure/plink/coil.SNP.fam" "$STRICT_PHENO" \
  "$PCA.eigenvec" "$ROOT/01_structure/summary" 5

KIN="$ROOT/01_structure/gemma/output/coil_SNP_kinship.cXX.txt"
if [[ ! -s "$KIN" ]]; then
  mkdir -p "$ROOT/01_structure/gemma"
  (
    cd "$ROOT/01_structure/gemma"
    "$GEMMA" -bfile "$LD" -gk 1 -o coil_SNP_kinship
  )
fi
test -s "$KIN"

COVAR="$ROOT/01_structure/summary/covariates.pc5.txt"
PHENO_VECTOR="$ROOT/01_structure/summary/phenotype.0_1.txt"

run_one() {
  local label=$1
  local result="$ROOT/02_snp_indel_gwas/gemma/output/${label}.coil.pc5.assoc.txt"
  if [[ ! -s "$result" ]]; then
    echo "[$(date '+%F %T')] GEMMA $label start"
    (
      cd "$ROOT/02_snp_indel_gwas/gemma"
      "$GEMMA" -bfile "$ROOT/01_structure/plink/coil.$label" \
        -k "$KIN" -c "$COVAR" -p "$PHENO_VECTOR" -lmm 4 -o "${label}.coil.pc5"
    )
  fi
  test -s "$result"
}

run_one SNP &
pid_snp=$!
run_one INDEL &
pid_indel=$!
wait "$pid_snp"
wait "$pid_indel"

"$PYTHON" "$ROOT/scripts/summarize_gwas.py" \
  --assoc "SNP=$ROOT/02_snp_indel_gwas/gemma/output/SNP.coil.pc5.assoc.txt" \
  --assoc "INDEL=$ROOT/02_snp_indel_gwas/gemma/output/INDEL.coil.pc5.assoc.txt" \
  --outdir "$ROOT/02_snp_indel_gwas/results" --prefix coil_SNP_INDEL

{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'samples\t125\n'
  printf 'phenotype\t0=non_spiral,1=spiral; weak_spiral_excluded\n'
  printf 'model\tGEMMA_LMM_Wald\n'
  printf 'kinship\tLD_pruned_SNP_centered\n'
  printf 'covariates\tintercept_plus_PC1-PC5\n'
  printf 'marker_qc\tMAF>=0.05; missingness<=0.20\n'
} > "$ROOT/summary/snp_indel_gwas.done"

echo "[$(date '+%F %T %Z')] strict coil SNP/INDEL GWAS completed"
