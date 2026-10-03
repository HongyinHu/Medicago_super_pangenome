#!/usr/bin/env bash
#SBATCH --job-name=gwas_chr23997_finalize
#SBATCH --partition=pLiu,pNormal
#SBATCH --qos=fast
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=path/to/project/N_4.pod_spiny/28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/logs/slurm_chr23997_finalize_%j.out
#SBATCH --error=path/to/project/N_4.pod_spiny/28_GWAS_spiny/08_Chr23997_targeted_genotyping_20260715/logs/slurm_chr23997_finalize_%j.err

set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/08_Chr23997_targeted_genotyping_20260715"
STRUCT="$ROOT/03_population_structure"
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
GEMMA=path/to/home/anaconda3/envs/biosofeware/bin/gemma

mkdir -p "$STAGE/results/vcf" "$STAGE/results/plink" "$STAGE/results/gemma/output" \
    "$STAGE/summary" "$STAGE/tmp"

completed=$(find "$STAGE/results/per_sample" -mindepth 2 -maxdepth 2 -name '*.done' -size +0c | wc -l)
if [[ "$completed" -ne 143 ]]; then
    echo "Expected 143 complete samples, found $completed" >&2
    exit 2
fi

"$PYTHON" "$STAGE/scripts/summarize_target_genotypes.py" \
    --manifest "$STAGE/inputs/target_samples.143.tsv" \
    --per-sample-dir "$STAGE/results/per_sample" \
    --summary-dir "$STAGE/summary" \
    --vcf-dir "$STAGE/results/vcf"

for mode in strict sensitivity; do
    input_vcf="$STAGE/results/vcf/Chr23997.target.${mode}.proxy.vcf"
    output_vcf="$input_vcf.gz"
    "$BCFTOOLS" view -Oz -o "$output_vcf" "$input_vcf"
    "$BCFTOOLS" index -f -t "$output_vcf"

    if [[ ! -s "$STAGE/summary/${mode}.gemma_eligible" ]]; then
        printf 'callset\t%s\nstatus\tnot_testable_no_polymorphic_called_site\n' "$mode" \
            > "$STAGE/summary/${mode}.gemma_status.tsv"
        continue
    fi

    prefix="$STAGE/results/plink/Chr23997.${mode}"
    "$PLINK" --vcf "$output_vcf" --double-id --chr-set 8 no-xy \
        --keep-allele-order --make-bed --out "$prefix" \
        >"$STAGE/logs/plink.${mode}.out" 2>"$STAGE/logs/plink.${mode}.err"
    awk 'BEGIN {OFS="\t"} {print $1, $2}' "$prefix.fam" > "$STAGE/tmp/${mode}.samples"
    awk 'BEGIN {OFS="\t"} {print $1, $2}' "$STRUCT/summary/analysis_samples.keep" \
        > "$STAGE/tmp/expected.samples"
    if ! cmp -s "$STAGE/tmp/${mode}.samples" "$STAGE/tmp/expected.samples"; then
        echo "Sample order mismatch for $mode" >&2
        diff -u "$STAGE/tmp/expected.samples" "$STAGE/tmp/${mode}.samples" >&2 || true
        exit 3
    fi
    [[ $(wc -l < "$prefix.bim") -eq 1 ]]

    (
        cd "$STAGE/results/gemma"
        "$GEMMA" -bfile "$prefix" \
            -k "$STRUCT/results/gemma/output/SNP_kinship.cXX.txt" \
            -c "$STRUCT/summary/covariates.pc5.txt" \
            -p "$STAGE/summary/${mode}.phenotype.called_only.txt" \
            -lmm 4 -o "Chr23997.${mode}.pc5"
    ) >"$STAGE/logs/gemma.${mode}.out" 2>"$STAGE/logs/gemma.${mode}.err"
    assoc="$STAGE/results/gemma/output/Chr23997.${mode}.pc5.assoc.txt"
    test -s "$assoc"
    {
        printf 'callset\t%s\n' "$mode"
        printf 'status\tcompleted\n'
        printf 'model\tGEMMA_LMM_Wald_complete_case\n'
        printf 'kinship\tSNP_LD_pruned_centered\n'
        printf 'covariates\tPC1-PC5\n'
        awk 'NR == 2 {printf "variant\t%s\nbeta\t%s\nse\t%s\np_wald\t%s\n", $2, $9, $10, $13}' "$assoc"
    } > "$STAGE/summary/${mode}.gemma_status.tsv"
done

{
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
    printf 'samples\t143\n'
    printf 'event\tChr4:88980831-88981052 DEL 221bp\n'
    printf 'main_callset\tstrict\n'
    printf 'association_primary\tGEMMA_LMM_complete_case_with_SNP_kinship_and_PC1-PC5_when_polymorphic\n'
    printf 'association_unadjusted\ttwo_sided_Fisher_exact\n'
    printf 'depth_policy\tQC_only_never_sufficient_alone_for_DEL\n'
    printf 'missing_policy\tinsufficient_or_conflicting_evidence_is_missing_not_0/0\n'
    printf 'graphTyper_note\tprevious_exact_target_output_had_zero_variant_records\n'
} > "$STAGE/summary/targeted_genotyping.done"

cat "$STAGE/summary/Chr23997.targeted_association.fisher.tsv"
for status in "$STAGE"/summary/*.gemma_status.tsv; do
    echo "--- $status"
    cat "$status"
done
