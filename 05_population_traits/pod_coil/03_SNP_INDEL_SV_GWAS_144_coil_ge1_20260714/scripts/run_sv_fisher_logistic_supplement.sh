#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
OUT="$ROOT/07_SV_fisher_logistic_supplement_20260715"
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
TABIX=path/to/home/anaconda3/envs/panpop/bin/tabix
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
export PATH=path/to/home/anaconda3/envs/panpop/bin:path/to/home/anaconda3/envs/biosofeware/bin:/usr/bin:/bin
export TMPDIR="$OUT/tmp"

mkdir -p "$OUT"/{filtered,model,qc,fisher,logistic,results,figures,logs,tmp,summary}
exec > >(tee -a "$OUT/logs/run_sv_fisher_logistic_supplement.log") 2>&1

# Start before the original MAF filter, then reproduce sample and marker missingness
# and MAF filtering in the retained 137-sample SV cohort.
PRE_QC="$ROOT/03_sv_gwas/plink/coil144.SV.GQ20.geno20"
COMMON_KEEP="$ROOT/01_structure/model/common_analysis.keep"
BASE="$OUT/filtered/SV.GQ20.mind20.maf005"
FINAL="$OUT/filtered/SV.GQ20.common.geno20.maf005"
RAW_VCF="$ROOT/03_sv_gwas/joint/coil144.SV.GQ20.vcf.gz"
[[ -s "$PRE_QC.bed" && -s "$RAW_VCF" ]] || { echo 'Missing established SV GWAS inputs' >&2; exit 2; }

if [[ ! -s "$BASE.bed" ]]; then
  "$PLINK" --bfile "$PRE_QC" --chr-set 8 no-xy --allow-no-sex \
    --mind 0.20 --maf 0.05 --make-bed --out "$BASE"
fi
if [[ ! -s "$FINAL.bed" ]]; then
  "$PLINK" --bfile "$BASE" --chr-set 8 no-xy --allow-no-sex --keep "$COMMON_KEEP" \
    --geno 0.20 --maf 0.05 --make-bed --out "$FINAL"
fi
"$PLINK" --bfile "$PRE_QC" --chr-set 8 no-xy --allow-no-sex --missing --out "$OUT/qc/SV.pre_final_qc"
"$PLINK" --bfile "$FINAL" --chr-set 8 no-xy --allow-no-sex --missing --out "$OUT/qc/SV.final_maf005_missing"

N=$(wc -l < "$FINAL.fam")
CASE=$(awk -F '\t' 'NR>1 && $2==1 {n++} END{print n+0}' "$ROOT/00_input/coil_ge1_phenotype.tsv")
(( N >= 100 )) || { echo "Too few SV samples after QC: $N" >&2; exit 3; }

"$PYTHON" "$ROOT/scripts/make_model_metadata.py" "$FINAL.fam" \
  "$ROOT/00_input/coil_ge1_phenotype.tsv" "$ROOT/01_structure/plink/coil144.SNP.pca.eigenvec" \
  "$OUT/model" SV 5

LOGISTIC="$OUT/logistic/SV.maf005.geno20.pc5"
if [[ ! -s "$LOGISTIC.assoc.logistic" ]]; then
  "$PLINK" --bfile "$FINAL" --chr-set 8 no-xy --allow-no-sex \
    --pheno "$OUT/model/SV.plink.pheno.tsv" --covar "$OUT/model/SV.plink.covar.tsv" \
    --covar-name PC1,PC2,PC3,PC4,PC5 --logistic hide-covar --ci 0.95 --out "$LOGISTIC"
fi

# Use raw GQ-filtered VCF genotypes to call ALT-event carriers, rather than the
# potentially minor-allele-reordered PLINK A1 coding used for logistic regression.
# SV IDs are not unique in the raw VCF, so filtering by PLINK marker ID would
# incorrectly retain multiple independent records. Reproduce the original
# two-stage PLINK QC at record level: --geno 0.20 across all 144 samples, then
# --geno 0.20 and MAF >=0.05 in the retained 137-sample cohort.
"$BCFTOOLS" query -l "$RAW_VCF" > "$OUT/qc/SV.all144.vcf_samples.txt"
[[ $(wc -l < "$OUT/qc/SV.all144.vcf_samples.txt") -eq 144 ]] || { echo 'Expected 144 raw VCF samples' >&2; exit 4; }
QUERY="$OUT/fisher/SV.twostage_qc.GT.tsv"
if [[ ! -s "$QUERY" ]]; then
  "$BCFTOOLS" view -m2 -M2 -Ou "$RAW_VCF" \
    | "$BCFTOOLS" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' > "$QUERY"
fi
BIM_N=$(wc -l < "$FINAL.bim")
SAMPLE_AUDIT=$(awk -F '\t' 'NR>1 {if($2==1)c++; else if($2==0)n++} END{printf "N=%d; case=%d; control=%d", c+n,c,n}' "$OUT/model/SV.model.tsv")
"$PYTHON" "$ROOT/scripts/analyze_sv_fisher_logistic.py" \
  --query "$QUERY" --all-samples "$OUT/qc/SV.all144.vcf_samples.txt" --phenotype "$OUT/model/SV.model.tsv" \
  --fisher-output "$OUT/fisher/SV.maf005.geno20.fisher.tsv" --logistic "$LOGISTIC.assoc.logistic" \
  --outdir "$OUT" --sample-audit "$SAMPLE_AUDIT" --maf 0.05 --max-initial-missing 0.20 --max-final-missing 0.20 --expected-records "$BIM_N"
FISHER_N=$(($(wc -l < "$OUT/fisher/SV.maf005.geno20.fisher.tsv") - 1))

{
  echo "completed	$(date '+%F %T %Z')"
  echo "input	GQ>=20 SV VCF"
  echo "stage1_marker_missingness_all144	<=0.20"
  echo "stage2_samples	137 common GWAS samples; marker missingness<=0.20"
  echo "MAF	>=0.05"
  echo "fisher	Two-sided ALT-carrier Fisher exact test"
  echo "logistic	PLINK additive logistic regression with PC1-PC5"
  echo "samples	$SAMPLE_AUDIT"
  echo "markers	$FISHER_N"
} > "$OUT/summary/SV_fisher_logistic.done"
