#!/usr/bin/env bash
set -euo pipefail

STAGE=path/to/project/N_4.pod_spiny/30_Chr23997_targeted_graphgenotype_GWAS_20260716
MANIFEST=$STAGE/inputs/target_samples.143.tsv
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
JOBS=${JOBS:-12}

mkdir -p "$STAGE/results/breakpoint_per_sample" "$STAGE/logs/breakpoint_all" "$STAGE/summary"
run_one() {
    local task="$1"
    local line task_index sample phenotype bam out
    line=$(awk -F '\t' -v task="$task" 'NR == task + 1 {print; exit}' "$MANIFEST")
    IFS=$'\t' read -r task_index sample phenotype bam <<< "$line"
    out="$STAGE/results/breakpoint_per_sample/${sample}.tsv"
    if [[ -s "$out" ]]; then
        return 0
    fi
    "$PYTHON" "$STAGE/scripts/chr23997_breakpoint_genotyper.py" \
        --bam "$bam" --sample "$sample" --out "$out"
}

for task in $(awk 'NR>1 {print $1}' "$MANIFEST"); do
    while (( $(jobs -pr | wc -l) >= JOBS )); do sleep 1; done
    run_one "$task" >"$STAGE/logs/breakpoint_all/task_${task}.out" 2>"$STAGE/logs/breakpoint_all/task_${task}.err" &
done
wait

"$PYTHON" "$STAGE/scripts/summarize_chr23997_breakpoint_calls.py" \
    --manifest "$MANIFEST" \
    --calls-dir "$STAGE/results/breakpoint_per_sample" \
    --out "$STAGE/summary/Chr23997_breakpoint_calls.143.tsv" \
    --summary "$STAGE/summary/Chr23997_breakpoint_fisher.143.tsv"
touch "$STAGE/summary/breakpoint_genotyping.done"
