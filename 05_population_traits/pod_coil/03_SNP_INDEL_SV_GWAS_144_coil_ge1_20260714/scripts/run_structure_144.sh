#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
PHENO=path/to/project/N_4.pod_spiny/21_resequence_phenotype_completed_20260708/results/resequence_sample_pod_traits_completed.tsv
AUDIT=path/to/project/N_4.pod_spiny/28_GWAS_spiny/01_input_audit/summary/phenotype_144.tsv
SNP_SOURCE=path/to/project/38.medicago_resequence/2.call_SNP_new/07_filter_snp/DP_6-85_miss_0.2.all_samples.biallelic.SNP.vcf.gz
INDEL_SOURCE=path/to/project/38.medicago_resequence/2.call_SNP_new/08_filter_indel/DP_6-85_miss_0.2.all_samples.biallelic.INDEL.vcf.gz
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python

mkdir -p "$ROOT"/{00_input,01_structure/{bcf,plink,model},02_snp_indel_gwas/{gmmat,logistic,results},03_sv_gwas,04_candidate_gene_peaks,logs,summary,tmp,software/Rlib}
LOG="$ROOT/logs/run_structure_144.log"
exec > >(tee -a "$LOG") 2>&1

echo "[$(date '+%F %T %Z')] 144-sample structure preparation start on $(hostname)"
"$PYTHON" "$ROOT/scripts/prepare_coil_ge1.py" "$PHENO" "$AUDIT" "$ROOT/00_input"

SAMPLES="$ROOT/00_input/samples_144.txt"
MAP="$ROOT/tmp/chrom_rename.tsv"
printf 'Chr1\t1\nChr2\t2\nChr3\t3\nChr4\t4\nChr5\t5\nChr6\t6\nChr7\t7\nChr8\t8\n' > "$MAP"

validate_source_samples() {
  local label=$1 source=$2
  "$BCFTOOLS" query -l "$source" | sort -u > "$ROOT/tmp/${label}.source.samples.sorted"
  sort -u "$SAMPLES" > "$ROOT/tmp/desired.samples.sorted"
  comm -23 "$ROOT/tmp/desired.samples.sorted" "$ROOT/tmp/${label}.source.samples.sorted" > "$ROOT/summary/${label}.missing_from_source.txt"
  [[ ! -s "$ROOT/summary/${label}.missing_from_source.txt" ]] || {
    echo "$label source VCF is missing requested samples" >&2
    cat "$ROOT/summary/${label}.missing_from_source.txt" >&2
    exit 2
  }
}

extract_and_qc() {
  local label=$1 source=$2
  local bcf="$ROOT/01_structure/bcf/coil144.${label}.bcf"
  local stage1="$ROOT/01_structure/plink/coil144.${label}.geno20"
  local final="$ROOT/01_structure/plink/coil144.${label}.twostage"
  validate_source_samples "$label" "$source"
  if [[ ! -s "$bcf" || ! -s "$bcf.csi" ]]; then
    "$BCFTOOLS" view -S "$SAMPLES" -r Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8 -Ou "$source" \
      | "$BCFTOOLS" annotate --rename-chrs "$MAP" --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ob -o "$bcf"
    "$BCFTOOLS" index -f "$bcf"
  fi
  [[ $("$BCFTOOLS" query -l "$bcf" | wc -l) -eq 144 ]]
  if [[ ! -s "$stage1.bed" ]]; then
    "$PLINK" --bcf "$bcf" --double-id --chr-set 8 no-xy --allow-extra-chr --geno 0.20 \
      --make-bed --out "$stage1"
  fi
  if [[ ! -s "$final.bed" ]]; then
    "$PLINK" --bfile "$stage1" --chr-set 8 no-xy --mind 0.20 --maf 0.05 \
      --make-bed --out "$final"
  fi
  local retained
  retained=$(wc -l < "$final.fam")
  (( retained >= 120 )) || { echo "$label retained only $retained samples" >&2; exit 3; }
}

extract_and_qc SNP "$SNP_SOURCE"
extract_and_qc INDEL "$INDEL_SOURCE"

for label in SNP INDEL; do
  awk '{print $2}' "$ROOT/01_structure/plink/coil144.${label}.twostage.fam" | sort -u \
    > "$ROOT/tmp/${label}.post_qc.ids"
done
comm -12 "$ROOT/tmp/SNP.post_qc.ids" "$ROOT/tmp/INDEL.post_qc.ids" \
  > "$ROOT/01_structure/model/common_analysis.ids"
COMMON_N=$(wc -l < "$ROOT/01_structure/model/common_analysis.ids")
(( COMMON_N >= 120 )) || { echo "Only $COMMON_N samples passed both SNP and INDEL QC" >&2; exit 4; }
awk 'NR==FNR {keep[$1]=1; next} ($2 in keep) {print $1, $2}' \
  "$ROOT/01_structure/model/common_analysis.ids" \
  "$ROOT/01_structure/plink/coil144.SNP.twostage.fam" \
  > "$ROOT/01_structure/model/common_analysis.keep"

printf 'class\tinput_samples\tclass_post_qc_samples\tcommon_analysis_samples\tpost_qc_markers\n' \
  > "$ROOT/summary/post_qc_sample_marker_counts.tmp.tsv"
for label in SNP INDEL; do
  source_prefix="$ROOT/01_structure/plink/coil144.${label}.twostage"
  analysis_prefix="$ROOT/01_structure/plink/coil144.${label}.analysis"
  if [[ ! -s "$analysis_prefix.bed" ]]; then
    "$PLINK" --bfile "$source_prefix" --chr-set 8 no-xy \
      --keep "$ROOT/01_structure/model/common_analysis.keep" \
      --make-bed --out "$analysis_prefix"
  fi
  [[ $(wc -l < "$analysis_prefix.fam") -eq $COMMON_N ]]
  printf '%s\t%s\t%s\t%s\t%s\n' "$label" 144 \
    "$(wc -l < "$source_prefix.fam")" "$COMMON_N" "$(wc -l < "$analysis_prefix.bim")" \
    >> "$ROOT/summary/post_qc_sample_marker_counts.tmp.tsv"
done
mv "$ROOT/summary/post_qc_sample_marker_counts.tmp.tsv" "$ROOT/summary/post_qc_sample_marker_counts.tsv"

SNP="$ROOT/01_structure/plink/coil144.SNP.analysis"
LD="$ROOT/01_structure/plink/coil144.SNP.ldpruned"
if [[ ! -s "$SNP.prune.in" ]]; then
  "$PLINK" --bfile "$SNP" --chr-set 8 no-xy --indep-pairwise 50 5 0.2 --out "$SNP"
fi
if [[ ! -s "$LD.bed" ]]; then
  "$PLINK" --bfile "$SNP" --extract "$SNP.prune.in" --chr-set 8 no-xy --make-bed --out "$LD"
fi

PCA="$ROOT/01_structure/plink/coil144.SNP.pca"
if [[ ! -s "$PCA.eigenvec" ]]; then
  "$PLINK" --bfile "$LD" --chr-set 8 no-xy --pca 10 --out "$PCA"
fi
GRM="$ROOT/01_structure/plink/coil144.SNP.grm"
if [[ ! -s "$GRM.rel" || ! -s "$GRM.rel.id" ]]; then
  "$PLINK" --bfile "$LD" --chr-set 8 no-xy --make-rel square --out "$GRM"
fi

for label in SNP INDEL; do
  "$PYTHON" "$ROOT/scripts/make_model_metadata.py" \
    "$ROOT/01_structure/plink/coil144.${label}.analysis.fam" \
    "$ROOT/00_input/coil_ge1_phenotype.tsv" "$PCA.eigenvec" \
    "$ROOT/01_structure/model" "$label" 5
  for chrom in $(seq 1 8); do
    chr_out="$ROOT/01_structure/plink/coil144.${label}.analysis.chr${chrom}"
    if [[ ! -s "$chr_out.bed" ]]; then
      "$PLINK" --bfile "$ROOT/01_structure/plink/coil144.${label}.analysis" \
        --chr-set 8 no-xy --chr "$chrom" --make-bed --out "$chr_out"
    fi
  done
done

bash "$ROOT/scripts/install_gmmat.sh"

{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'input_samples\t144\n'
  printf 'phenotype\taverage_coil>=1 spiral; average_coil<1 non-spiral\n'
  printf 'composition\t103 spiral; 41 non-spiral\n'
  printf 'common_analysis_samples\t%s\n' "$COMMON_N"
  printf 'marker_qc\tstage1 missingness<=0.20; stage2 sample missingness<=0.20 and MAF>=0.05\n'
  printf 'structure\tLD-pruned SNP GRM and PC1-PC5\n'
} > "$ROOT/summary/structure_144.done"
echo "[$(date '+%F %T %Z')] structure preparation complete"
