#!/usr/bin/env bash
set -euo pipefail
RUN=path/to/project/N_4.pod_spiny/20_Chr23997_SNP_INDEL_SV_peak_20260708
GWAS=path/to/project/N_4.pod_spiny/08_GWAS_pod_spine_SNP_SV_20260706
SV=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706
BCF=path/to/home/anaconda3/envs/panpop/bin/bcftools
TABIX=path/to/home/anaconda3/envs/panpop/bin/tabix
PY=path/to/home/anaconda3/envs/panpop/bin/python
PLINK2=path/to/home/anaconda3/envs/call_snp/bin/plink2
PHENO=$GWAS/input/pod_spine_plink2_header.pheno
COVAR=$GWAS/qc/sativa182_Msa_pod_spine_snp.maf005.pca.eigenvec
INDEL=path/to/data
DELLY=$SV/results/delly_DEL_cohort/phenotyped99.DEL.genotyped.vcf.gz
mkdir -p "$RUN/results" "$RUN/summary" "$RUN/logs" "$RUN/figures"
exec > >(tee -a "$RUN/logs/run_chr23997_peak_assoc.log") 2>&1
date; hostname
if [ ! -s "$RUN/results/indel.Chr4_88_90Mb.biallelic.vcf.gz" ]; then
  $BCF view -r Chr4:88000000-90000000 -m2 -M2 -Oz -o "$RUN/results/indel.Chr4_88_90Mb.biallelic.vcf.gz" "$INDEL"
  $TABIX -f -p vcf "$RUN/results/indel.Chr4_88_90Mb.biallelic.vcf.gz"
fi
if [ ! -s "$RUN/results/indel.Chr4_88_90Mb.pgen" ]; then
  $PLINK2 --vcf "$RUN/results/indel.Chr4_88_90Mb.biallelic.vcf.gz" --double-id --allow-extra-chr --vcf-half-call missing --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 --make-pgen --maf 0.05 --out "$RUN/results/indel.Chr4_88_90Mb"
fi
$PLINK2 --pfile "$RUN/results/indel.Chr4_88_90Mb" --allow-extra-chr --pheno "$PHENO" --pheno-name pod_spine --glm hide-covar --out "$RUN/results/indel.Chr4_88_90Mb.no_covar"
$PLINK2 --pfile "$RUN/results/indel.Chr4_88_90Mb" --allow-extra-chr --pheno "$PHENO" --pheno-name pod_spine --covar "$COVAR" --covar-name PC1,PC2,PC3,PC4,PC5 --glm hide-covar --out "$RUN/results/indel.Chr4_88_90Mb.PC5"
if [ ! -s "$RUN/results/dellyDEL.Chr4_88_90Mb.vcf.gz" ]; then
  $BCF view -r Chr4:88000000-90000000 -Oz -o "$RUN/results/dellyDEL.Chr4_88_90Mb.vcf.gz" "$DELLY"
  $TABIX -f -p vcf "$RUN/results/dellyDEL.Chr4_88_90Mb.vcf.gz"
fi
if [ ! -s "$RUN/results/dellyDEL.Chr4_88_90Mb.pgen" ]; then
  $PLINK2 --vcf "$RUN/results/dellyDEL.Chr4_88_90Mb.vcf.gz" --double-id --allow-extra-chr --vcf-half-call missing --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 --make-pgen --maf 0.05 --out "$RUN/results/dellyDEL.Chr4_88_90Mb"
fi
$PLINK2 --pfile "$RUN/results/dellyDEL.Chr4_88_90Mb" --allow-extra-chr --pheno "$PHENO" --pheno-name pod_spine --glm hide-covar --out "$RUN/results/dellyDEL.Chr4_88_90Mb.no_covar"
$PLINK2 --pfile "$RUN/results/dellyDEL.Chr4_88_90Mb" --allow-extra-chr --pheno "$PHENO" --pheno-name pod_spine --covar "$COVAR" --covar-name PC1,PC2,PC3,PC4,PC5 --glm hide-covar --out "$RUN/results/dellyDEL.Chr4_88_90Mb.PC5"
$PY "$RUN/scripts/summarize_gwas_region.py" \
  SNP_no_covar=logistic="$GWAS/results/sativa182_pod_spine_snp_gwas.corrected.no_covar.pod_spine.glm.logistic" \
  SNP_PC5=logistic="$GWAS/results/sativa182_pod_spine_snp_gwas.corrected.PC5.pod_spine.glm.logistic" \
  INDEL_no_covar=logistic="$RUN/results/indel.Chr4_88_90Mb.no_covar.pod_spine.glm.logistic" \
  INDEL_PC5=logistic="$RUN/results/indel.Chr4_88_90Mb.PC5.pod_spine.glm.logistic" \
  DELLY_DEL_no_covar=logistic="$RUN/results/dellyDEL.Chr4_88_90Mb.no_covar.pod_spine.glm.logistic" \
  DELLY_DEL_PC5=logistic="$RUN/results/dellyDEL.Chr4_88_90Mb.PC5.pod_spine.glm.logistic" \
  > "$RUN/summary/Chr23997_region_SNP_INDEL_SV_top_hits.tsv"
CAND=$SV/results/Chr23997_candidate_genotype/Chr23997_DEL_221.depth_ratio.valid_phenotyped.with_gt.tsv
if [ -s "$CAND" ]; then
  $PY "$RUN/scripts/test_chr23997_candidate_del.py" "$CAND" "$RUN/summary/Chr23997_DEL_221_direct_association.txt" || true
fi
date
touch "$RUN/summary/chr23997_peak_assoc.done"
