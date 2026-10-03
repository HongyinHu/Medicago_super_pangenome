#!/usr/bin/env bash
set -euo pipefail

LABEL=${1:?Usage: run_modelA_maf_marker.sh SNP_or_INDEL maf chromosome}
MAF=${2:?Usage: run_modelA_maf_marker.sh SNP_or_INDEL maf chromosome}
CHR=${SLURM_ARRAY_TASK_ID:-${3:-}}
[[ "$LABEL" == SNP || "$LABEL" == INDEL ]] || { echo "LABEL must be SNP or INDEL" >&2; exit 2; }
[[ "$CHR" =~ ^[1-8]$ ]] || { echo "Chromosome must be 1-8" >&2; exit 2; }
awk -v maf="$MAF" 'BEGIN { exit !(maf > 0 && maf < 0.5) }'
MAF_TAG=$(awk -v maf="$MAF" 'BEGIN { printf "maf%03d", int(maf * 100 + 0.5) }')

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
OUT=$ROOT/06_maf_sensitivity_modelA_20260715
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript
export GWAS_ROOT="$ROOT"
export R_LIBS_USER="$ROOT/software/Rlib"
export TMPDIR="$OUT/tmp/$MAF_TAG/${LABEL}_chr${CHR}"

mkdir -p "$TMPDIR" "$OUT/$MAF_TAG/bfiles/$LABEL" "$OUT/$MAF_TAG/raw/$LABEL" "$OUT/logs" "$OUT/summary/tasks"
LOG="$OUT/logs/${MAF_TAG}.${LABEL}.chr${CHR}.log"
exec > >(tee -a "$LOG") 2>&1

test -s "$ROOT/summary/structure_144.done"
SOURCE="$ROOT/01_structure/plink/coil144.${LABEL}.analysis.chr${CHR}"
MODEL="$ROOT/01_structure/model/${LABEL}.model.tsv"
GRM="$ROOT/01_structure/plink/coil144.SNP.grm"
FILTERED="$OUT/$MAF_TAG/bfiles/$LABEL/${LABEL}.chr${CHR}.${MAF_TAG}"
RESULT="$OUT/$MAF_TAG/raw/$LABEL/${LABEL}.chr${CHR}.Model_A_${MAF_TAG}.gmmat.score.tsv"
DONE="$OUT/summary/tasks/${MAF_TAG}.${LABEL}.chr${CHR}.done"
test -s "$SOURCE.bed"

echo "[$(date '+%F %T %Z')] $LABEL Chr${CHR}, Model A, MAF >= $MAF on $(hostname)"
if [[ ! -s "$FILTERED.bed" || ! -s "$FILTERED.bim" || ! -s "$FILTERED.fam" ]]; then
  "$PLINK" --bfile "$SOURCE" --chr-set 8 no-xy --allow-no-sex --maf "$MAF" --make-bed --out "$FILTERED"
fi
test -s "$FILTERED.bed"
test -s "$FILTERED.bim"
test -s "$FILTERED.fam"
MARKERS=$(wc -l < "$FILTERED.bim")
if [[ "$MARKERS" -eq 0 ]]; then
  echo "No markers remain after MAF filtering" >&2
  exit 1
fi

if [[ ! -s "$RESULT" ]]; then
  "$RSCRIPT" --vanilla "$ROOT/scripts/run_gmmat_gwas_maf.R" \
    "$FILTERED" "$MODEL" "$GRM.rel" "$GRM.rel.id" "$RESULT" \
    "${LABEL}_Model_A_${MAF_TAG}_Chr${CHR}" 5 "${SLURM_CPUS_PER_TASK:-2}" "$MAF"
fi
test -s "$RESULT"
awk -v min_maf="$MAF" 'NR > 1 { af=$8 + 0; maf=(af <= 0.5 ? af : 1-af); if (maf + 1e-9 < min_maf) { print "Below requested MAF:", $2, maf > "/dev/stderr"; exit 1 } }' "$RESULT"

tmp_done="$DONE.tmp.$$"
{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'model\tA: GMMAT binomial mixed model; PC1-PC5; SNP GRM\n'
  printf 'variant_class\t%s\n' "$LABEL"
  printf 'chromosome\t%s\n' "$CHR"
  printf 'maf_cutoff\t%s\n' "$MAF"
  printf 'samples\t%s\n' "$(wc -l < "$FILTERED.fam")"
  printf 'markers_after_plink_maf_filter\t%s\n' "$MARKERS"
  printf 'score_file\t%s\n' "$RESULT"
} > "$tmp_done"
mv "$tmp_done" "$DONE"
echo "[$(date '+%F %T %Z')] complete"
