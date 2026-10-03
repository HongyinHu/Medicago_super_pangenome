#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
mkdir -p "$ROOT/02_snp_indel_gwas/results" "$ROOT/logs" "$ROOT/summary"
exec > >(tee -a "$ROOT/logs/run_snp_indel_finalize.log") 2>&1
for label in SNP INDEL; do
  for chrom in $(seq 1 8); do
    test -s "$ROOT/summary/${label}.chr${chrom}.done"
  done
  combined="$ROOT/02_snp_indel_gwas/gmmat/${label}.coil_ge1.pc5.gmmat.score.tsv"
  files=()
  for chrom in $(seq 1 8); do
    files+=("$ROOT/02_snp_indel_gwas/gmmat/${label}.chr${chrom}.coil_ge1.pc5.gmmat.score.tsv")
  done
  awk 'FNR==1 && NR!=1 {next} {print}' "${files[@]}" > "$combined"
  test -s "$combined"
  {
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
    printf 'variant_class\t%s\n' "$label"
    printf 'chromosome_tasks\t8\n'
    printf 'parallelization\tone single-core PLINK scan per chromosome\n'
  } > "$ROOT/summary/${label}_gwas.done"
done
"$PYTHON" "$ROOT/scripts/summarize_gmmat.py" \
  --assoc "SNP=$ROOT/02_snp_indel_gwas/gmmat/SNP.coil_ge1.pc5.gmmat.score.tsv" \
  --assoc "INDEL=$ROOT/02_snp_indel_gwas/gmmat/INDEL.coil_ge1.pc5.gmmat.score.tsv" \
  --outdir "$ROOT/02_snp_indel_gwas/results" --prefix coil144_SNP_INDEL
printf 'completed\t%s\n' "$(date '+%F %T %Z')" > "$ROOT/summary/snp_indel_finalize.done"
