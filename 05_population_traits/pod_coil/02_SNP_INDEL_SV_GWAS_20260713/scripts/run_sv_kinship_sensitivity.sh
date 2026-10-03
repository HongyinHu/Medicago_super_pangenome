#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/02_SNP_INDEL_SV_GWAS_20260713
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
SV_PLINK="$ROOT/03_sv_gwas/plink/coil.SV.GQ20.twostage"
OUTROOT="$ROOT/03_sv_gwas/kinship_only"
PREFIX=SV.coil.kinship_only.gq20.twostage

mkdir -p "$OUTROOT/gemma" "$OUTROOT/results"
[[ $(wc -l < "$SV_PLINK.fam") -eq 119 ]]

if [[ ! -s "$OUTROOT/gemma/output/$PREFIX.assoc.txt" ]]; then
  (
    cd "$OUTROOT/gemma"
    "$GEMMA" -bfile "$SV_PLINK" \
      -k "$ROOT/03_sv_gwas/summary/SV_kinship.reordered.cXX.txt" \
      -p "$ROOT/03_sv_gwas/summary/SV.phenotype.0_1.txt" \
      -lmm 4 -o "$PREFIX"
  )
fi

"$PYTHON" "$ROOT/scripts/summarize_gwas.py" \
  --assoc "SV=$OUTROOT/gemma/output/$PREFIX.assoc.txt" \
  --outdir "$OUTROOT/results" --prefix coil_SV_kinship_only

printf 'completed\t%s\nsamples\t119\nmodel\tGEMMA_LMM_Wald; SNP_kinship; no_PCs\n' \
  "$(date '+%F %T %Z')" > "$OUTROOT/kinship_only.done"
