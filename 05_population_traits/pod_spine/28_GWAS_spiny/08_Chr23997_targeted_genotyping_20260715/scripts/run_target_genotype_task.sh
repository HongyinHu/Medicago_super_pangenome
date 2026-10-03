#!/usr/bin/env bash
#SBATCH --job-name=gwas_chr23997_gt
#SBATCH --partition=pLiu,pNormal
#SBATCH --qos=fast
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=1-143%96
#SBATCH --output=path/to/project/N_4.pod_spiny/28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/logs/slurm_chr23997_gt_%A_%a.out
#SBATCH --error=path/to/project/N_4.pod_spiny/28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/logs/slurm_chr23997_gt_%A_%a.err

set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/08_Chr23997_targeted_genotyping_20260715"
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
SMOOVE=path/to/home/anaconda3/envs/smoove/bin/smoove
DELLY=path/to/home/anaconda3/envs/delly/bin/delly
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python

export PATH=path/to/home/anaconda3/envs/smoove/bin:path/to/home/anaconda3/envs/delly/bin:path/to/home/anaconda3/envs/panpop/bin:/usr/bin:/bin

manifest="$STAGE/inputs/target_samples.143.tsv"
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
line=$(awk -F '\t' -v target="$task_id" 'NR == target + 1 {print; exit}' "$manifest")
if [[ -z "$line" ]]; then
    echo "No manifest row for task $task_id" >&2
    exit 2
fi
IFS=$'\t' read -r manifest_index sample phenotype bam <<< "$line"
if [[ "$manifest_index" != "$task_id" ]]; then
    echo "Manifest/task mismatch: $manifest_index != $task_id" >&2
    exit 2
fi

outdir="$STAGE/results/per_sample/$sample"
logdir="$STAGE/logs/per_sample/$sample"
mkdir -p "$outdir" "$logdir"
done_file="$outdir/$sample.done"

validate_outputs() {
    test -s "$outdir/${sample}-smoove.genotyped.vcf.gz"
    test -s "$outdir/${sample}-smoove.genotyped.vcf.gz.csi"
    [[ $("$BCFTOOLS" view -H "$outdir/${sample}-smoove.genotyped.vcf.gz" | wc -l) -eq 1 ]]
    test -s "$outdir/${sample}.delly.bcf"
    test -s "$outdir/${sample}.delly.bcf.csi"
    "$BCFTOOLS" view -h "$outdir/${sample}.delly.bcf" >/dev/null
    test -s "$outdir/${sample}.evidence.tsv"
    [[ $(wc -l < "$outdir/${sample}.evidence.tsv") -eq 2 ]]
}

if [[ -s "$done_file" ]] && validate_outputs; then
    echo "[$(date '+%F %T')] $sample already complete"
    exit 0
fi

test -s "$bam"
if [[ ! -s "$bam.bai" && ! -s "${bam%.bam}.bai" ]]; then
    echo "Missing BAM index for $bam" >&2
    exit 4
fi

echo "[$(date '+%F %T')] sample=$sample phenotype=$phenotype node=$(hostname) start"
"$SMOOVE" genotype --name "$sample" --outdir "$outdir" --fasta "$REF" \
    --removepr --processes 1 \
    --vcf "$STAGE/inputs/Chr23997.intron2_DEL.svtyper.vcf.gz" "$bam" \
    >"$logdir/smoove.out" 2>"$logdir/smoove.err"

"$DELLY" call -g "$REF" -v "$STAGE/inputs/Chr23997.intron2_DEL.delly.bcf" \
    -o "$outdir/${sample}.delly.bcf" "$bam" \
    >"$logdir/delly.out" 2>"$logdir/delly.err"
"$BCFTOOLS" index -f "$outdir/${sample}.delly.bcf"

"$PYTHON" "$STAGE/scripts/extract_chr23997_bam_evidence.py" \
    --bam "$bam" --sample "$sample" --out "$outdir/${sample}.evidence.tsv" \
    >"$logdir/evidence.out" 2>"$logdir/evidence.err"

validate_outputs
{
    printf 'sample\t%s\n' "$sample"
    printf 'phenotype\t%s\n' "$phenotype"
    printf 'node\t%s\n' "$(hostname)"
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
} > "$done_file"
echo "[$(date '+%F %T')] sample=$sample complete"
