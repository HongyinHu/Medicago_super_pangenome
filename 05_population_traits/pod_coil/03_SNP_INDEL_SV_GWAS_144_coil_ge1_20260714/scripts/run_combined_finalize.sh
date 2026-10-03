#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
mkdir -p "$ROOT/04_candidate_gene_peaks" "$ROOT/logs" "$ROOT/summary"
exec > >(tee -a "$ROOT/logs/run_combined_finalize.log") 2>&1
for label in SNP INDEL SV; do test -s "$ROOT/summary/${label}_gwas.done"; done
"$PYTHON" "$ROOT/scripts/summarize_gmmat.py" \
  --assoc "SNP=$ROOT/02_snp_indel_gwas/gmmat/SNP.coil_ge1.pc5.gmmat.score.tsv" \
  --assoc "INDEL=$ROOT/02_snp_indel_gwas/gmmat/INDEL.coil_ge1.pc5.gmmat.score.tsv" \
  --assoc "SV=$ROOT/03_sv_gwas/gmmat/SV.coil_ge1.pc5.gmmat.score.tsv" \
  --outdir "$ROOT/04_candidate_gene_peaks" --prefix coil144_SNP_INDEL_SV
printf 'completed\t%s\n' "$(date '+%F %T %Z')" > "$ROOT/summary/pipeline.done"

