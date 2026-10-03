#!/usr/bin/env bash
set -euo pipefail

RUN=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STAGE="$RUN/05_sv_genotyping"
SOURCE=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706/medicago144_3caller_20260709
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
KEEP="$RUN/03_population_structure/summary/analysis_samples.keep"
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
SMOOVE=path/to/home/anaconda3/envs/smoove/bin/smoove
PYTHON=/usr/bin/python3
FILTER_SCRIPT="$STAGE/scripts/filter_smoove_sites.py"
PAIR_FILTER_SCRIPT="$STAGE/scripts/remove_orphan_bnd.py"
PAIR_VALIDATOR="$STAGE/scripts/validate_bnd_pairs.py"
FILTER_JOBS=${FILTER_JOBS:-8}
MIN_DISCOVERY_SAMPLES=${MIN_DISCOVERY_SAMPLES:-2}
MAX_SVLEN=${MAX_SVLEN:-1000000}
PRIMARY_CONTIGS=Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8

export PATH=path/to/home/anaconda3/envs/smoove/bin:path/to/home/anaconda3/envs/panpop/bin:/usr/bin:/bin
mkdir -p "$STAGE"/{candidate_sites,filtered_discovery,logs,summary,tmp}
LOG="$STAGE/logs/run_sv_site_discovery.log"
exec > >(tee -a "$LOG") 2>&1

echo "[$(date '+%F %T %Z')] SV site discovery start on $(hostname)"
[[ -s "$KEEP" && -s "$REF" && -s "$REF.fai" && -s "$FILTER_SCRIPT" ]]
[[ -s "$PAIR_FILTER_SCRIPT" && -s "$PAIR_VALIDATOR" ]]

MANIFEST="$STAGE/summary/sv_analysis_manifest.tsv"
VCF_LIST="$STAGE/tmp/filtered_discovery_vcfs.list"
printf 'sample_id\tbam\tdiscovery_vcf\tfiltered_vcf\n' > "$MANIFEST"
: > "$VCF_LIST"

while read -r fid sid rest; do
    [[ -n "$sid" ]] || continue
    bam=$(awk -F '\t' -v s="$sid" 'NR>1 && $1==s {print $6; exit}' "$SOURCE/input/medicago144_sample_bams.tsv")
    src="$SOURCE/results/smoove_per_sample/$sid/$sid-smoove.genotyped.vcf.gz"
    out="$STAGE/filtered_discovery/$sid.smoove.discovery.filtered.vcf.gz"
    [[ -s "$bam" && ( -s "$bam.bai" || -s "${bam%.bam}.bai" ) ]] || { echo "Missing BAM/index for $sid" >&2; exit 1; }
    [[ -s "$src" && ( -s "$src.csi" || -s "$src.tbi" ) ]] || { echo "Missing Smoove VCF/index for $sid" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\n' "$sid" "$bam" "$src" "$out" >> "$MANIFEST"
    printf '%s\n' "$out" >> "$VCF_LIST"
done < "$KEEP"

samples=$(( $(wc -l < "$MANIFEST") - 1 ))
[[ "$samples" -eq 143 ]] || { echo "Expected 143 samples, found $samples" >&2; exit 1; }

filter_one() {
    local sid=$1 src=$2 out=$3
    if [[ -s "$out" && ( -s "$out.csi" || -s "$out.tbi" ) ]] && "$BCFTOOLS" view -h "$out" >/dev/null 2>&1; then
        return 0
    fi
    local tmp="${out%.gz}.tmp.vcf.gz"
    rm -f "$tmp" "$tmp.csi" "$tmp.tbi"
    # PRPOS/PREND are intentionally retained here because svtools lmerge
    # requires LUMPY breakpoint-probability vectors. They are removed only
    # after official Smoove/SVtools clustering has completed.
    "$BCFTOOLS" view -r "$PRIMARY_CONTIGS" \
        -i 'INFO/SU>=5 && FORMAT/GQ>=20 && QUAL>=20 && (INFO/SVTYPE="BND" || (abs(INFO/SVLEN)>=50 && abs(INFO/SVLEN)<=1000000))' \
        "$src" -Oz -o "$tmp"
    "$BCFTOOLS" index -f "$tmp"
    mv "$tmp" "$out"
    if [[ -s "$tmp.csi" ]]; then mv "$tmp.csi" "$out.csi"; fi
    if [[ -s "$tmp.tbi" ]]; then mv "$tmp.tbi" "$out.tbi"; fi
}
export -f filter_one
export BCFTOOLS PRIMARY_CONTIGS

running=0
while IFS=$'\t' read -r sid bam src out; do
    [[ "$sid" == sample_id ]] && continue
    filter_one "$sid" "$src" "$out" &
    running=$((running + 1))
    if (( running >= FILTER_JOBS )); then
        wait -n
        running=$((running - 1))
    fi
done < "$MANIFEST"
wait

COUNTS="$STAGE/summary/filtered_discovery_counts.tsv"
printf 'sample_id\trecords\n' > "$COUNTS"
while IFS=$'\t' read -r sid bam src out; do
    [[ "$sid" == sample_id ]] && continue
    n=$("$BCFTOOLS" index -n "$out")
    printf '%s\t%s\n' "$sid" "$n" >> "$COUNTS"
done < "$MANIFEST"

RAW="$STAGE/candidate_sites/medicago143.sites.raw.vcf.gz"
if [[ ! -s "$RAW" ]]; then
    rm -f "$STAGE/candidate_sites/medicago143.lsort.vcf" \
          "$STAGE/candidate_sites/medicago143.sites.vcf.gz" \
          "$STAGE/candidate_sites/medicago143.sites.vcf.gz.csi" \
          "$STAGE/candidate_sites/medicago143.sites.vcf.gz.tbi"
    mapfile -t vcfs < "$VCF_LIST"
    "$SMOOVE" merge --name medicago143 --outdir "$STAGE/candidate_sites" --fasta "$REF" "${vcfs[@]}"
    [[ -s "$STAGE/candidate_sites/medicago143.sites.vcf.gz" ]]
    mv "$STAGE/candidate_sites/medicago143.sites.vcf.gz" "$RAW"
    rm -f "$STAGE/candidate_sites/medicago143.lsort.vcf"
fi

PRELIM="$STAGE/tmp/medicago143.sites.filtered.preliminary.vcf"
PLAIN="$STAGE/tmp/medicago143.sites.filtered.paired.vcf"
SMOOVE_SITES="$STAGE/candidate_sites/medicago143.sites.filtered.nsamp${MIN_DISCOVERY_SAMPLES}.smoove.vcf.gz"
FINAL="$STAGE/candidate_sites/medicago143.sites.filtered.nsamp${MIN_DISCOVERY_SAMPLES}.vcf.gz"
FILTER_SUMMARY="$STAGE/summary/candidate_site_filter_counts.tsv"
PAIR_SUMMARY="$STAGE/summary/orphan_bnd_filter_counts.tsv"
PAIR_VALIDATION="$STAGE/summary/bnd_pair_validation.smoove_order.tsv"
"$PYTHON" "$FILTER_SCRIPT" --input "$RAW" --output "$PRELIM" \
    --min-discovery-samples "$MIN_DISCOVERY_SAMPLES" --max-svlen "$MAX_SVLEN" --summary "$FILTER_SUMMARY"
"$PYTHON" "$PAIR_FILTER_SCRIPT" --input "$PRELIM" --output "$PLAIN" --summary "$PAIR_SUMMARY"
rm -f "$SMOOVE_SITES"
bgzip -@ 4 -c "$PLAIN" > "$SMOOVE_SITES"
"$PYTHON" "$PAIR_VALIDATOR" "$SMOOVE_SITES" --output "$PAIR_VALIDATION" --require-adjacent
# svtools lmerge can emit non-contiguous chromosome blocks, especially when
# interchromosomal BND records are retained. Sort after site filtering so the
# indexed QC VCF is deterministic. Smoove genotype uses SMOOVE_SITES instead,
# where BND mates remain adjacent in the original lmerge order.
rm -f "$FINAL" "$FINAL.tbi" "$FINAL.csi"
"$BCFTOOLS" sort -m 2G -T "$STAGE/tmp/bcftools-sort.XXXXXX" \
    -Oz -o "$FINAL" "$PLAIN"
tabix -f -p vcf "$FINAL"
"$BCFTOOLS" view -h "$FINAL" >/dev/null
records=$("$BCFTOOLS" index -n "$FINAL")
rm -f "$PRELIM" "$PLAIN"

cat > "$STAGE/summary/sv_site_discovery.done" <<EOF
completed\t$(date '+%F %T %Z')
samples\t$samples
candidate_sites\t$records
min_single_sample_SU\t5
min_single_sample_GQ\t20
min_discovery_samples\t$MIN_DISCOVERY_SAMPLES
primary_contigs\tChr1-Chr8
interval_svlen\t50-$MAX_SVLEN
types\tDEL,DUP,INV,BND
candidate_vcf\t$FINAL
smoove_genotype_vcf\t$SMOOVE_SITES
EOF
echo "[$(date '+%F %T %Z')] completed with $records candidate sites"
