#!/usr/bin/env bash
set -euo pipefail
ROOT=path/to/project
RUN="$ROOT/N_5.pod_coil/01_SHP_homology_and_coil_cosegregation"
DATA="$ROOT/N_5.pod_coil/00_data"
BCF=path/to/home/anaconda3/envs/panpop/bin/bcftools
VCF="$DATA/SV_SNP/DP_6-85_miss_0.2.all_samples.biallelic.SNP.vcf.gz"
PHENO="$DATA/144_phenotype/medicago144_pod_spiral_spine_traits.tsv"
"$BCF" query -l "$VCF" > "$RUN/population/vcf.samples"
for spec in 'Chr14120 Chr3:53194386-53205316' 'Chr24229 Chr4:91171307-91180034'; do
  set -- $spec
  gene=$1
  region=$2
  "$BCF" query -r "$region" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$VCF" > "$RUN/population/${gene}.SNP.genotypes.tsv"
  python3 "$RUN/scripts/association.py" population --query "$RUN/population/${gene}.SNP.genotypes.tsv" \
    --samples "$RUN/population/vcf.samples" --phenotype "$PHENO" --gene "$gene" \
    --output "$RUN/tables/${gene}.population_SNP_coil_association.tsv"
done
touch "$RUN/population_stats.done"
