#!/usr/bin/env bash
set -euo pipefail

RUN=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$RUN/05_sv_genotyping"
MANIFEST="$STAGE/summary/sv_analysis_manifest.tsv"
SITES_INDEXED="$STAGE/candidate_sites/medicago143.sites.filtered.nsamp2.vcf.gz"
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
SUMMARY="$STAGE/summary/sv_joint_genotype_validation.tsv"
COMPLETE="$STAGE/summary/sv_joint_genotypes.validated.done"
FAILED="$STAGE/summary/sv_joint_genotypes.validation.failed.tsv"
TMP="$STAGE/tmp/sv_joint_genotype_validation.$$"

mkdir -p "$STAGE/tmp" "$STAGE/summary"
rm -f "$COMPLETE" "$FAILED"
trap 'rm -f "$TMP"' EXIT

expected_records=$("$BCFTOOLS" index -n "$SITES_INDEXED")
expected_samples=$(awk 'END {print NR - 1}' "$MANIFEST")
printf 'sample_id\tvcf\tsize_bytes\tgzip_ok\theader_ok\tindex_ok\trecords\tsample_name_ok\tstatus\n' > "$TMP"

failures=0
checked=0
while IFS=$'\t' read -r sid bam discovery_vcf filtered_vcf; do
    [[ "$sid" == "sample_id" ]] && continue
    checked=$((checked + 1))
    out="$STAGE/joint_genotype_per_sample/$sid/$sid-smoove.genotyped.vcf.gz"
    done_file="$STAGE/joint_genotype_per_sample/$sid/$sid.done"
    size=0
    gzip_ok=0
    header_ok=0
    index_ok=0
    records=NA
    sample_ok=0
    status=FAIL

    if [[ -s "$out" ]]; then
        size=$(stat -c '%s' "$out")
        gzip -t "$out" >/dev/null 2>&1 && gzip_ok=1
        "$BCFTOOLS" view -h "$out" >/dev/null 2>&1 && header_ok=1
        if [[ -s "$out.csi" ]]; then
            records=$("$BCFTOOLS" index -n "$out" 2>/dev/null || printf NA)
            [[ "$records" == "$expected_records" ]] && index_ok=1
        fi
        sample_name=$("$BCFTOOLS" query -l "$out" 2>/dev/null || true)
        [[ "$sample_name" == "$sid" ]] && sample_ok=1
    fi

    if [[ "$gzip_ok" == 1 && "$header_ok" == 1 && "$index_ok" == 1 && "$sample_ok" == 1 ]]; then
        status=PASS
        printf 'completed\t%s\nhost\t%s\nsample\t%s\nrecords\t%s\nmode\trecovered_valid_genotype_no_duphold\n' \
            "$(date '+%F %T %Z')" "$(hostname)" "$sid" "$records" > "$done_file"
    else
        failures=$((failures + 1))
        rm -f "$done_file"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sid" "$out" "$size" "$gzip_ok" "$header_ok" "$index_ok" "$records" "$sample_ok" "$status" >> "$TMP"
done < "$MANIFEST"

mv "$TMP" "$SUMMARY"
trap - EXIT

if [[ "$checked" != "$expected_samples" || "$failures" != 0 ]]; then
    awk -F '\t' 'NR == 1 || $NF != "PASS"' "$SUMMARY" > "$FAILED"
    printf 'Validation failed: checked=%s expected=%s failures=%s\n' "$checked" "$expected_samples" "$failures" >&2
    exit 1
fi

printf 'completed\t%s\nhost\t%s\nsamples\t%s\nrecords_per_sample\t%s\nmode\tjoint_genotype_without_duphold\n' \
    "$(date '+%F %T %Z')" "$(hostname)" "$checked" "$expected_records" > "$COMPLETE"
printf 'Validation passed: samples=%s records_per_sample=%s\n' "$checked" "$expected_records"
