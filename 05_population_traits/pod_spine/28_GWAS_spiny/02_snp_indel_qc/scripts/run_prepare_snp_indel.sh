#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$ROOT/02_snp_indel_qc"
AUDIT="$ROOT/01_input_audit/summary"
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
TABIX=path/to/home/anaconda3/envs/panpop/bin/tabix
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
PYTHON=/usr/bin/python3
SNP_IN=path/to/project/38.medicago_resequence/2.call_SNP_new/07_filter_snp/DP_6-85_miss_0.2.all_samples.biallelic.SNP.vcf.gz
INDEL_IN=path/to/project/38.medicago_resequence/2.call_SNP_new/08_filter_indel/DP_6-85_miss_0.2.all_samples.biallelic.INDEL.vcf.gz
SAMPLES="$AUDIT/sample_order.txt"
PHENO="$AUDIT/phenotype_144.tsv"

mkdir -p "$STAGE/results/vcf" "$STAGE/results/plink" "$STAGE/logs" "$STAGE/summary" "$STAGE/tmp"

MAP="$STAGE/chrom_rename.tsv"
printf 'Chr1\t1\nChr2\t2\nChr3\t3\nChr4\t4\nChr5\t5\nChr6\t6\nChr7\t7\nChr8\t8\n' > "$MAP"

prepare_one() {
  local label=$1
  local input_vcf=$2
  local output_vcf="$STAGE/results/vcf/144.${label}.Chr1-8.vcf.gz"
  local prefix="$STAGE/results/plink/144.${label}.common"
  local done="$STAGE/summary/${label}.prepare.done"

  if [[ ! -s "$output_vcf" || ! -s "$output_vcf.tbi" ]]; then
    echo "[$(date '+%F %T')] Subsetting $label to the 144 phenotype samples and Chr1-Chr8"
    "$BCFTOOLS" view \
      --regions Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8 \
      --samples-file "$SAMPLES" --force-samples \
      --output-type u "$input_vcf" \
      | "$BCFTOOLS" annotate \
          --rename-chrs "$MAP" \
          --set-id '%CHROM:%POS:%REF:%FIRST_ALT' \
          --output-type z --output "$output_vcf"
    "$TABIX" -f -p vcf "$output_vcf"
  fi

  if [[ ! -s "$prefix.bed" || ! -s "$prefix.bim" || ! -s "$prefix.fam" ]]; then
    echo "[$(date '+%F %T')] Converting $label to PLINK and applying common-marker QC"
    "$PLINK" \
      --vcf "$output_vcf" \
      --double-id \
      --chr-set 8 no-xy \
      --mind 0.20 \
      --geno 0.20 \
      --maf 0.05 \
      --make-bed \
      --out "$prefix"
    "$PYTHON" "$STAGE/scripts/set_fam_phenotype.py" "$PHENO" "$prefix.fam"
  fi

  "$PLINK" --bfile "$prefix" --chr-set 8 no-xy --freq --missing --out "$STAGE/summary/144.${label}.common"
  {
    printf 'label\t%s\n' "$label"
    printf 'samples\t%s\n' "$(wc -l < "$prefix.fam")"
    printf 'variants\t%s\n' "$(wc -l < "$prefix.bim")"
    printf 'maf_threshold\t0.05\n'
    printf 'variant_missingness_threshold\t0.20\n'
    printf 'sample_missingness_threshold\t0.20\n'
    printf 'HWE_filter\tnot_applied_multi_species_panel\n'
    printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  } > "$done"
}

prepare_one SNP "$SNP_IN"
prepare_one INDEL "$INDEL_IN"

touch "$STAGE/summary/prepare_snp_indel.done"
echo "[$(date '+%F %T')] SNP/INDEL preparation completed"
