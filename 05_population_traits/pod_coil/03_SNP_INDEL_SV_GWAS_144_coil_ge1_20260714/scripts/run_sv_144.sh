#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
SOURCE=path/to/project/N_4.pod_spiny/28_GWAS_spiny/05_sv_genotyping
PLINK=path/to/home/anaconda3/envs/biosofeware/bin/plink
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
TABIX=path/to/home/anaconda3/envs/panpop/bin/tabix
SMOOVE=path/to/home/anaconda3/envs/smoove/bin/smoove
PYTHON=path/to/home/anaconda3/envs/panpop/bin/python
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript
export PATH=path/to/home/anaconda3/envs/smoove/bin:path/to/home/anaconda3/envs/panpop/bin:path/to/home/anaconda3/envs/biosofeware/bin:/usr/bin:/bin
export GWAS_ROOT="$ROOT"
export R_LIBS_USER="$ROOT/software/Rlib"
export TMPDIR="$ROOT/tmp/R_SV"

mkdir -p "$ROOT/03_sv_gwas"/{joint,plink,model,gmmat,logistic,results} "$ROOT/logs" "$ROOT/summary" "$ROOT/tmp" "$TMPDIR"
LOG="$ROOT/logs/run_sv_144.log"
exec > >(tee -a "$LOG") 2>&1
test -s "$ROOT/summary/structure_144.done"
test -s "$ROOT/summary/qh3_sv_genotype.done"

VCF_LIST="$ROOT/tmp/sv_genotype_144.list"
awk -F '\t' 'NR>1 {print $1}' "$SOURCE/summary/sv_analysis_manifest.tsv" | while read -r sid; do
  vcf="$SOURCE/joint_genotype_per_sample/$sid/$sid-smoove.genotyped.vcf.gz"
  [[ -s "$vcf" ]] || { echo "Missing existing SV genotype: $vcf" >&2; exit 2; }
  printf '%s\n' "$vcf"
done > "$VCF_LIST"
printf '%s\n' "$ROOT/03_sv_gwas/qh3_genotype/qh3-smoove.genotyped.vcf.gz" >> "$VCF_LIST"
[[ $(wc -l < "$VCF_LIST") -eq 144 ]]

JOINT="$ROOT/03_sv_gwas/joint/medicago144.smoove.square.vcf.gz"
if [[ ! -s "$JOINT" ]]; then
  mapfile -t vcfs < "$VCF_LIST"
  "$SMOOVE" paste --name medicago144 --outdir "$ROOT/03_sv_gwas/joint" "${vcfs[@]}"
fi
if [[ ! -s "$JOINT.csi" && ! -s "$JOINT.tbi" ]]; then
  "$BCFTOOLS" index -f "$JOINT"
fi
[[ $("$BCFTOOLS" query -l "$JOINT" | wc -l) -eq 144 ]]
"$BCFTOOLS" query -l "$JOINT" | sort -u > "$ROOT/tmp/SV.joint.samples.sorted"
sort -u "$ROOT/00_input/samples_144.txt" > "$ROOT/tmp/desired.samples.sorted"
comm -3 "$ROOT/tmp/desired.samples.sorted" "$ROOT/tmp/SV.joint.samples.sorted" > "$ROOT/summary/SV_joint_sample_difference.txt"
[[ ! -s "$ROOT/summary/SV_joint_sample_difference.txt" ]]

MAP="$ROOT/tmp/chrom_rename.tsv"
printf 'Chr1\t1\nChr2\t2\nChr3\t3\nChr4\t4\nChr5\t5\nChr6\t6\nChr7\t7\nChr8\t8\n' > "$MAP"
VCF="$ROOT/03_sv_gwas/joint/coil144.SV.GQ20.vcf.gz"
if [[ ! -s "$VCF" || ! -s "$VCF.tbi" ]]; then
  "$BCFTOOLS" view -S "$ROOT/00_input/samples_144.txt" -r Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8 -Ou "$JOINT" \
    | "$BCFTOOLS" filter -S . -e 'FMT/GQ<20' -Ou \
    | "$BCFTOOLS" annotate --rename-chrs "$MAP" --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o "$VCF"
  "$TABIX" -f -p vcf "$VCF"
fi

STAGE1="$ROOT/03_sv_gwas/plink/coil144.SV.GQ20.geno20"
FINAL_QC="$ROOT/03_sv_gwas/plink/coil144.SV.GQ20.twostage"
FINAL="$ROOT/03_sv_gwas/plink/coil144.SV.GQ20.analysis"
if [[ ! -s "$STAGE1.bed" ]]; then
  "$PLINK" --vcf "$VCF" --double-id --chr-set 8 no-xy --allow-extra-chr --biallelic-only strict \
    --geno 0.20 --make-bed --out "$STAGE1"
fi
if [[ ! -s "$FINAL_QC.bed" ]]; then
  "$PLINK" --bfile "$STAGE1" --chr-set 8 no-xy --mind 0.20 --maf 0.05 \
    --make-bed --out "$FINAL_QC"
fi
if [[ ! -s "$FINAL.bed" ]]; then
  "$PLINK" --bfile "$FINAL_QC" --chr-set 8 no-xy \
    --keep "$ROOT/01_structure/model/common_analysis.keep" \
    --make-bed --out "$FINAL"
fi
SV_N=$(wc -l < "$FINAL.fam")
(( SV_N >= 100 )) || { echo "SV QC retained only $SV_N samples" >&2; exit 3; }

PCA="$ROOT/01_structure/plink/coil144.SNP.pca.eigenvec"
"$PYTHON" "$ROOT/scripts/make_model_metadata.py" "$FINAL.fam" \
  "$ROOT/00_input/coil_ge1_phenotype.tsv" "$PCA" "$ROOT/03_sv_gwas/model" SV 5

OUT="$ROOT/03_sv_gwas/gmmat/SV.coil_ge1.pc5.gmmat.score.tsv"
GRM="$ROOT/01_structure/plink/coil144.SNP.grm"
if [[ ! -s "$OUT" ]]; then
  "$RSCRIPT" --vanilla "$ROOT/scripts/run_gmmat_gwas.R" \
    "$FINAL" "$ROOT/03_sv_gwas/model/SV.model.tsv" "$GRM.rel" "$GRM.rel.id" \
    "$OUT" SV 5 "${SLURM_CPUS_PER_TASK:-8}"
fi
test -s "$OUT"

SENS="$ROOT/03_sv_gwas/logistic/SV.coil_ge1.pc5"
if [[ ! -s "$SENS.assoc.logistic" ]]; then
  "$PLINK" --bfile "$FINAL" --chr-set 8 no-xy --allow-no-sex \
    --pheno "$ROOT/03_sv_gwas/model/SV.plink.pheno.tsv" \
    --covar "$ROOT/03_sv_gwas/model/SV.plink.covar.tsv" \
    --covar-name PC1,PC2,PC3,PC4,PC5 --logistic hide-covar --ci 0.95 --out "$SENS"
fi

"$PYTHON" "$ROOT/scripts/summarize_gmmat.py" --assoc "SV=$OUT" \
  --outdir "$ROOT/03_sv_gwas/results" --prefix coil144_SV
{
  printf 'completed\t%s\n' "$(date '+%F %T %Z')"
  printf 'joint_input_samples\t144\n'
  printf 'post_qc_samples\t%s\n' "$SV_N"
  printf 'post_qc_markers\t%s\n' "$(wc -l < "$FINAL.bim")"
  printf 'genotype_qc\tGQ>=20; marker missingness<=0.20; sample missingness<=0.20; MAF>=0.05\n'
  printf 'primary_model\tGMMAT binomial logistic mixed model; SNP GRM; PC1-PC5\n'
} > "$ROOT/summary/SV_gwas.done"
