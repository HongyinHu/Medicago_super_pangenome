#!/usr/bin/env bash
set -euo pipefail

STAGE=path/to/project/N_4.pod_spiny/30_Chr23997_targeted_graphgenotype_GWAS_20260716
REFERENCE=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
BWA=path/to/home/anaconda3/envs/biosofeware/bin/bwa-mem2

"$PYTHON" "$STAGE/scripts/prepare_chr23997_haplotypes.py" \
    --reference "$REFERENCE" \
    --out "$STAGE/inputs/Chr23997_REF_DEL221.haplotypes.fa"
"$BWA" index "$STAGE/inputs/Chr23997_REF_DEL221.haplotypes.fa"

# Four phenotype/state strata from the prior caller are represented:
# phenotype 0 old REF, phenotype 0 old DEL, phenotype 1 old REF,
# and phenotype 1 old DEL.
for task in 4 7 14 26 31 40 49 63 64 92 114 115 116 120 126 140; do
    bash "$STAGE/scripts/run_chr23997_haplotype_task.sh" "$task"
done

"$PYTHON" "$STAGE/scripts/summarize_chr23997_haplotype_calls.py" \
    --manifest "$STAGE/inputs/target_samples.143.tsv" \
    --per-sample-dir "$STAGE/results/per_sample" \
    --out "$STAGE/summary/pilot_genotypes.tsv"
