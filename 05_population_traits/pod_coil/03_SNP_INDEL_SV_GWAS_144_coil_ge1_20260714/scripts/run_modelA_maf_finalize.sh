#!/usr/bin/env bash
set -euo pipefail

MAF=${1:?Usage: run_modelA_maf_finalize.sh maf}
MAF_TAG=$(awk -v maf="$MAF" 'BEGIN { printf "maf%03d", int(maf * 100 + 0.5) }')
ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
OUT=$ROOT/06_maf_sensitivity_modelA_20260715
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
RESULTS=$OUT/$MAF_TAG/results
RAW=$OUT/$MAF_TAG/raw
mkdir -p "$RESULTS" "$OUT/$MAF_TAG/combined" "$OUT/tmp/matplotlib" "$OUT/logs"
export MPLCONFIGDIR="$OUT/tmp/matplotlib"
LOG="$OUT/logs/${MAF_TAG}.finalize.log"
exec > >(tee -a "$LOG") 2>&1

for label in SNP INDEL; do
  for chr in $(seq 1 8); do
    test -s "$OUT/summary/tasks/${MAF_TAG}.${label}.chr${chr}.done"
  done
  first="$RAW/$label/${label}.chr1.Model_A_${MAF_TAG}.gmmat.score.tsv"
  combined="$OUT/$MAF_TAG/combined/${label}.Model_A_${MAF_TAG}.gmmat.score.tsv"
  head -n 1 "$first" > "$combined"
  for chr in $(seq 1 8); do
    tail -n +2 "$RAW/$label/${label}.chr${chr}.Model_A_${MAF_TAG}.gmmat.score.tsv" >> "$combined"
  done
  test -s "$combined"
done
test -s "$OUT/summary/tasks/${MAF_TAG}.SV.done"
SV_COMBINED="$OUT/$MAF_TAG/combined/SV.Model_A_${MAF_TAG}.gmmat.score.tsv"
cp -f "$RAW/SV/SV.Model_A_${MAF_TAG}.gmmat.score.tsv" "$SV_COMBINED"

PREFIX=Model_A_MAF${MAF_TAG#maf}
"$PYTHON" "$ROOT/scripts/summarize_maf_cutoff.py" \
  --assoc "SNP=$OUT/$MAF_TAG/combined/SNP.Model_A_${MAF_TAG}.gmmat.score.tsv" \
  --assoc "INDEL=$OUT/$MAF_TAG/combined/INDEL.Model_A_${MAF_TAG}.gmmat.score.tsv" \
  --assoc "SV=$SV_COMBINED" \
  --outdir "$RESULTS" --prefix "$PREFIX" --maf "$MAF"
for required in "$RESULTS/$PREFIX.variant_summary.tsv" "$RESULTS/$PREFIX.qq.png" "$RESULTS/$PREFIX.manhattan.png" "$RESULTS/$PREFIX.assoc.tsv.gz"; do
  test -s "$required"
done
printf 'completed\t%s\nmaf_cutoff\t%s\nprefix\t%s\n' "$(date '+%F %T %Z')" "$MAF" "$PREFIX" > "$OUT/summary/${MAF_TAG}.finalized.done"
