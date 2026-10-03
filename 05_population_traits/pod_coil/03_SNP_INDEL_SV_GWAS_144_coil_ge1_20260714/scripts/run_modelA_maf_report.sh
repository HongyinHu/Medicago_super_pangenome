#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
OUT=$ROOT/06_maf_sensitivity_modelA_20260715
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
mkdir -p "$OUT/results" "$OUT/model_comparison/freq" "$OUT/logs" "$OUT/tmp/matplotlib"
export MPLCONFIGDIR="$OUT/tmp/matplotlib"
LOG="$OUT/logs/report.log"
exec > >(tee -a "$LOG") 2>&1

for tag in maf005 maf010 maf015; do
  test -s "$OUT/summary/${tag}.finalized.done"
done
test -s "$OUT/summary/diagnostics.done"

"$PYTHON" "$ROOT/scripts/compare_maf_sensitivity.py" \
  --cutoff "0.05=$OUT/maf005/results" \
  --cutoff "0.10=$OUT/maf010/results" \
  --cutoff "0.15=$OUT/maf015/results" \
  --outdir "$OUT/results"

for label in SNP INDEL; do
  source="$ROOT/01_structure/plink/coil144.${label}.analysis"
  "$PLINK" --bfile "$source" --chr-set 8 no-xy --allow-no-sex --freq --out "$OUT/model_comparison/freq/$label"
done
"$PLINK" --bfile "$ROOT/03_sv_gwas/plink/coil144.SV.GQ20.analysis" --chr-set 8 no-xy --allow-no-sex --freq --out "$OUT/model_comparison/freq/SV"
"$PYTHON" "$ROOT/scripts/compare_gwas_models.py" --root "$ROOT" --freqdir "$OUT/model_comparison/freq" --outdir "$OUT/model_comparison"
"$PYTHON" "$ROOT/scripts/build_gwas_diagnostic_report.py" --outdir "$OUT"

for required in \
  "$OUT/results/Model_A_MAF_sensitivity.tsv" \
  "$OUT/results/Model_A_MAF_sensitivity_lambda.png" \
  "$OUT/results/GWAS_diagnostic_report.md" \
  "$OUT/model_comparison/Model_ABC_top20_overlap.tsv" \
  "$OUT/model_comparison/Model_ABC_pvalue_correlation.tsv"; do
  test -s "$required"
done
printf 'completed\t%s\n' "$(date '+%F %T %Z')" > "$OUT/summary/modelA_maf_sensitivity.done"
