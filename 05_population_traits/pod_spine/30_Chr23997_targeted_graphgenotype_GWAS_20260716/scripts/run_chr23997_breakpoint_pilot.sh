#!/usr/bin/env bash
set -euo pipefail

STAGE=path/to/project/N_4.pod_spiny/30_Chr23997_targeted_graphgenotype_GWAS_20260716
MANIFEST=$STAGE/inputs/target_samples.143.tsv
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
TASKS=(4 7 14 26 31 40 49 63 64 92 114 115 116 120 126 140)

mkdir -p "$STAGE/results/breakpoint_per_sample" "$STAGE/summary"
for task in "${TASKS[@]}"; do
    line=$(awk -F '\t' -v task="$task" 'NR == task + 1 {print; exit}' "$MANIFEST")
    IFS=$'\t' read -r task_index sample phenotype bam <<< "$line"
    out="$STAGE/results/breakpoint_per_sample/${sample}.tsv"
    "$PYTHON" "$STAGE/scripts/chr23997_breakpoint_genotyper.py" \
        --bam "$bam" --sample "$sample" --out "$out"
    printf '%s\t%s\t%s\n' "$task" "$sample" "$phenotype"
done | tee "$STAGE/logs/run_chr23997_breakpoint_pilot.log"

{
    head -1 "$STAGE/results/breakpoint_per_sample/$(awk -F '\t' 'NR==5 {print $2}' "$MANIFEST").tsv"
    for task in "${TASKS[@]}"; do
        sample=$(awk -F '\t' -v task="$task" 'NR == task + 1 {print $2; exit}' "$MANIFEST")
        tail -n +2 "$STAGE/results/breakpoint_per_sample/${sample}.tsv"
    done
} > "$STAGE/summary/Chr23997_breakpoint_pilot.tsv"
