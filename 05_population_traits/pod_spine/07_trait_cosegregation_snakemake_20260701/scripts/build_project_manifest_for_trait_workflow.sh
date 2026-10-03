#!/usr/bin/env bash
set -euo pipefail

PKG=path/to/project/N_4.pod_spiny/07_trait_cosegregation_snakemake_20260701
SVROOT=path/to/project/N_3.call_SV/01.read_based_dualref_hifi

mkdir -p "$PKG/config"

vcf_out="$PKG/config/project_medicago_vcf_manifest.tsv"
bam_out="$PKG/config/project_medicago_bam_manifest.tsv"

printf "sample\tref\tcaller\tvcf\n" > "$vcf_out"
if [ -d "$SVROOT/03_per_sample" ]; then
  find "$SVROOT/03_per_sample" -type f -name "*.vcf.gz" | sort | while read -r vcf; do
    rel=${vcf#"$SVROOT/03_per_sample/"}
    ref=${rel%%/*}
    rest=${rel#*/}
    caller=${rest%%/*}
    base=$(basename "$vcf")
    sample=${base%%.$caller.vcf.gz}
    if [ "$sample" = "$base" ]; then
      sample=${base%%.vcf.gz}
    fi
    printf "%s\t%s\t%s\t%s\n" "$sample" "$ref" "$caller" "$vcf"
  done >> "$vcf_out"
fi

printf "sample\tref\tbam\n" > "$bam_out"
find "$SVROOT" -type f -name "*.bam" \
  ! -path "*/tools/*" ! -path "*/sofeware/*" ! -path "*/example/*" 2>/dev/null | sort | while read -r bam; do
    base=$(basename "$bam")
    sample=${base%%.bam}
    ref=unknown
    case "$bam" in
      *"/Msa/"*) ref=Msa ;;
      *"/R108/"*) ref=R108 ;;
      *".Msa."*) ref=Msa ;;
      *".R108."*) ref=R108 ;;
    esac
    printf "%s\t%s\t%s\n" "$sample" "$ref" "$bam"
  done >> "$bam_out"

echo "WROTE $vcf_out"
wc -l "$vcf_out"
echo "WROTE $bam_out"
wc -l "$bam_out"
