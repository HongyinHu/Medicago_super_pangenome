#!/usr/bin/env bash
set -euo pipefail

RUN=path/to/project/N_5.pod_coil/08_SHP_multiallelic_candidate_20260717
SOURCE=path/to/project/38.medicago_resequence/2.call_SNP_new/07_filter_snp/DP_6-85_miss_0.2.all_samples.SNP.vcf.gz
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
BCFTOOLS=path/to/software/samtools/bcftools-1.17/bcftools
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript
GWAS=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
THREADS=${THREADS:-8}
MODE=${1:-prepare}
REGIONS=Chr3:53194386-53205316,Chr4:91171307-91180034
MULTI=$RUN/vcf/coil144.SHP_AG.plusminus2kb.multiallelic.SNP.vcf.gz
SPLIT=$RUN/vcf/coil144.SHP_AG.plusminus2kb.split.SNP.vcf.gz
BFILE=$RUN/plink/coil143.SHP_AG.split
MODEL=$GWAS/01_structure/model/SNP.model.tsv
GRM=$GWAS/01_structure/plink/coil144.SNP.grm.rel
GRM_ID=$GWAS/01_structure/plink/coil144.SNP.grm.rel.id
KEEP=$GWAS/01_structure/model/common_analysis.keep

mkdir -p "$RUN"/{inputs,vcf,results,logs,summary,gmmat,plink,scripts,tests}
exec > >(tee -a "$RUN/logs/${MODE}.log") 2>&1

[[ -s "$SOURCE" && -s "$SOURCE.tbi" ]]
[[ -s "$REF" && -s "$REF.fai" ]]
[[ -x "$BCFTOOLS" ]]
[[ -x "$PLINK" ]]
[[ -x "$RSCRIPT" ]]

prepare() {
    echo "[$(date '+%F %T %Z')] extracting non-biallelic candidate SNPs on $(hostname)"
    "$BCFTOOLS" view --threads "$THREADS" -r "$REGIONS" -v snps -Oz -o "$MULTI" "$SOURCE"
    "$BCFTOOLS" index --threads "$THREADS" -t "$MULTI"
    "$BCFTOOLS" norm --threads "$THREADS" -f "$REF" -m-any -Ou "$MULTI" | \
        "$BCFTOOLS" view --threads "$THREADS" -v snps -e 'ALT="*"' -Oz -o "$SPLIT"
    "$BCFTOOLS" index --threads "$THREADS" -t "$SPLIT"

    samples=$("$BCFTOOLS" query -l "$SPLIT" | wc -l)
    records=$("$BCFTOOLS" index -n "$SPLIT")
    [[ "$samples" -eq 184 ]]
    [[ "$records" -gt 0 ]]

    python3 "$RUN/scripts/candidate_variant_tools.py" audit \
        --vcf "$SPLIT" \
        --phenotype "$GWAS/00_input/coil_ge1_phenotype.tsv" \
        --analysis-keep "$GWAS/01_structure/model/common_analysis.keep" \
        --out-genotypes "$RUN/results/functional_allele_genotypes.tsv" \
        --out-summary "$RUN/results/functional_allele_summary.tsv" \
        --expected-vcf-samples 184

    printf 'samples\t%s\nrecords\t%s\nsource_vcf\t%s\nsplit_vcf\t%s\n' \
        "$samples" "$records" "$SOURCE" "$SPLIT" > "$RUN/summary/candidate_vcf_audit.tsv"
    touch "$RUN/summary/candidate_vcf_audit.done"
    echo "[$(date '+%F %T %Z')] candidate VCF preparation complete"
}

association() {
    [[ -s "$SPLIT" && -s "$SPLIT.tbi" ]]
    echo "[$(date '+%F %T %Z')] making low-frequency candidate PLINK input on $(hostname)"
    "$PLINK" \
        --vcf "$SPLIT" \
        --double-id \
        --allow-extra-chr \
        --keep "$KEEP" \
        --indiv-sort file "$KEEP" \
        --geno 0.20 \
        --maf 0.005 \
        --set-missing-var-ids '@:#:$1:$2' \
        --new-id-max-allele-len 100 \
        --make-bed \
        --out "$BFILE"
    [[ $(wc -l < "$BFILE.fam") -eq 143 ]]
    [[ -s "$BFILE.bed" && -s "$BFILE.bim" ]]

    GWAS_ROOT="$GWAS" "$RSCRIPT" "$RUN/scripts/run_candidate_gmmat.R" \
        "$BFILE" "$MODEL" "$GRM" "$GRM_ID" \
        "$RUN/gmmat/candidate_pc1_5.tsv" candidate_pc1_5 5 "$THREADS" 0.005
    GWAS_ROOT="$GWAS" "$RSCRIPT" "$RUN/scripts/run_candidate_gmmat.R" \
        "$BFILE" "$MODEL" "$GRM" "$GRM_ID" \
        "$RUN/gmmat/candidate_pc1_2.tsv" candidate_pc1_2 2 "$THREADS" 0.005

    printf 'samples\t%s\nvariants\t%s\nmaf_min\t0.005\n' \
        "$(wc -l < "$BFILE.fam")" "$(wc -l < "$BFILE.bim")" > "$RUN/summary/candidate_plink_audit.tsv"
    echo "[$(date '+%F %T %Z')] candidate GMMAT association complete"
}

report() {
    [[ -s "$RUN/results/functional_allele_summary.tsv" ]]
    [[ -s "$BFILE.bim" ]]
    [[ -s "$RUN/gmmat/candidate_pc1_5.tsv" ]]
    [[ -s "$RUN/gmmat/candidate_pc1_2.tsv" ]]
    python3 "$RUN/scripts/build_candidate_report.py" \
        --functional-summary "$RUN/results/functional_allele_summary.tsv" \
        --bim "$BFILE.bim" \
        --pc1-5 "$RUN/gmmat/candidate_pc1_5.tsv" \
        --pc1-2 "$RUN/gmmat/candidate_pc1_2.tsv" \
        --out-tsv "$RUN/results/functional_allele_association_summary.tsv" \
        --out-md "$RUN/summary/candidate_association_report.md"
    echo "[$(date '+%F %T %Z')] candidate association report complete"
}

case "$MODE" in
    prepare) prepare ;;
    association) association ;;
    report) report ;;
    *) echo "Usage: $0 {prepare|association|report}" >&2; exit 2 ;;
esac
