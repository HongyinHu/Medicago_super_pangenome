#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
cd "$BASE"
mkdir -p results/delly_DEL_cohort summary logs
expected=$(wc -l < input/valid_phenotyped_bams.for_delly.tsv)
done_count=$(find results/delly_DEL_genotype_per_sample -name '*.delly.DEL.geno.done' 2>/dev/null | wc -l)
bcf_count=$(find results/delly_DEL_genotype_per_sample -name '*.delly.DEL.geno.bcf' 2>/dev/null | wc -l)
echo "expected=$expected done=$done_count bcf=$bcf_count"
if [[ "$done_count" -ne "$expected" || "$bcf_count" -ne "$expected" ]]; then
  echo "not all genotype outputs complete" >&2
  exit 3
fi
find results/delly_DEL_genotype_per_sample -name '*.delly.DEL.geno.bcf' | sort -V > input/delly_DEL_genotyped_bcf.list
outbcf=results/delly_DEL_cohort/phenotyped99.DEL.genotyped.bcf
outvcf=results/delly_DEL_cohort/phenotyped99.DEL.genotyped.vcf.gz
rm -f "$outbcf" "$outbcf.csi" "$outvcf" "$outvcf.tbi"
echo "[$(date '+%F %T')] bcftools merge start"
"$BCFTOOLS" merge -m id -Ob -o "$outbcf" -l input/delly_DEL_genotyped_bcf.list
"$BCFTOOLS" index -f "$outbcf"
"$BCFTOOLS" view -Oz -o "$outvcf" "$outbcf"
"$BCFTOOLS" index -t -f "$outvcf"
{
  echo -e "expected_samples\t$expected"
  echo -e "genotyped_bcf\t$bcf_count"
  echo -e "sites_count\t$($BCFTOOLS view -H results/delly_DEL_cohort/phenotyped99.DEL.sites.bcf | wc -l)"
  echo -e "merged_records\t$($BCFTOOLS view -H "$outvcf" | wc -l)"
  echo -e "sample_count\t$($BCFTOOLS query -l "$outvcf" | wc -l)"
  echo -e "out_vcf\t$BASE/$outvcf"
  echo -e "time\t$(date '+%F %T %Z')"
} > summary/delly_DEL_genotyped.summary.tsv
cat summary/delly_DEL_genotyped.summary.tsv
