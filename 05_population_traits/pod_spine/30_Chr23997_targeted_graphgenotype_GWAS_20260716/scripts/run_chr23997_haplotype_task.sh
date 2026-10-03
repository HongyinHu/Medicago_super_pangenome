#!/usr/bin/env bash
#SBATCH --job-name=chr23997_hapgt
#SBATCH --partition=pLiu,pNormal
#SBATCH --qos=fast
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail

STAGE=path/to/project/N_4.pod_spiny/30_Chr23997_targeted_graphgenotype_GWAS_20260716
REFERENCE=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
BWA=path/to/home/anaconda3/envs/biosofeware/bin/bwa-mem2
MANIFEST=${MANIFEST:-$STAGE/inputs/target_samples.143.tsv}
TASK_ID=${1:-${SLURM_ARRAY_TASK_ID:?provide task index or SLURM_ARRAY_TASK_ID}}
PIPELINE_VERSION=2

line=$(awk -F '\t' -v task="$TASK_ID" 'NR == task + 1 {print; exit}' "$MANIFEST")
if [[ -z "$line" ]]; then
    echo "No manifest row for task $TASK_ID" >&2
    exit 2
fi
IFS=$'\t' read -r task_index sample phenotype bam <<< "$line"
[[ "$task_index" == "$TASK_ID" ]] || { echo "Manifest/task mismatch" >&2; exit 2; }
[[ -s "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 3; }

out="$STAGE/results/per_sample/${sample}.tsv"
done_file="$STAGE/results/per_sample/${sample}.done"
mkdir -p "$STAGE/results/per_sample" "$STAGE/tmp/$sample" "$STAGE/logs"

if [[ -s "$done_file" && -s "$out" ]] && grep -Fqx $'pipeline_version\t'"$PIPELINE_VERSION" "$done_file"; then
    echo "$(date '+%F %T') ${sample} already complete"
    exit 0
fi

"$PYTHON" "$STAGE/scripts/chr23997_haplotype_genotyper.py" \
    --bam "$bam" \
    --reference "$REFERENCE" \
    --haplotypes "$STAGE/inputs/Chr23997_REF_DEL221.haplotypes.fa" \
    --bwa "$BWA" \
    --sample "$sample" \
    --out "$out" \
    --tmpdir "$STAGE/tmp/$sample"

{
    printf 'sample\t%s\n' "$sample"
    printf 'phenotype\t%s\n' "$phenotype"
    printf 'task_index\t%s\n' "$TASK_ID"
    printf 'pipeline_version\t%s\n' "$PIPELINE_VERSION"
    printf 'node\t%s\n' "$(hostname)"
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
} > "$done_file"
