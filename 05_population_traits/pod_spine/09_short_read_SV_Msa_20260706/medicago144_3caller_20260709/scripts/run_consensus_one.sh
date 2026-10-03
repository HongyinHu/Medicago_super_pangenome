#!/usr/bin/env bash
set -euo pipefail
MED=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706/medicago144_3caller_20260709
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
BGZIP=path/to/home/anaconda3/envs/panpop/bin/bgzip
TABIX=path/to/home/anaconda3/envs/panpop/bin/tabix
JAVA=path/to/home/anaconda3/envs/panpop/bin/java
JASMINE=path/to/project/N_3.call_SV/01.read_based_dualref_hifi/05_dualref_panSV_20260626/tools/jasmine.jar
TASK=${SLURM_ARRAY_TASK_ID:-${1:-}}
[ -n "$TASK" ] || { echo "Need task id" >&2; exit 2; }
LINE=$(awk -F'\t' -v n="$TASK" 'NR==n+1{print}' "$MED/input/medicago144_sample_bams.tsv")
[ -n "$LINE" ] || { echo "No line for task $TASK" >&2; exit 2; }
SID=$(echo "$LINE" | cut -f1)
NORMDIR="$MED/results/normalized/$SID"
OUTDIR="$MED/results/consensus_per_sample/$SID"
TMPDIR="$MED/tmp/consensus_$SID"
mkdir -p "$NORMDIR" "$OUTDIR" "$TMPDIR"
DONE="$OUTDIR/$SID.consensus_2of3.done"
OUTVCF="$OUTDIR/$SID.shortread_3caller.consensus2.vcf.gz"
if [ -s "$DONE" ] && [ -s "$OUTVCF" ] && [ -s "$OUTVCF.tbi" ]; then echo "skip consensus $SID"; exit 0; fi
rm -f "$DONE" "$OUTVCF" "$OUTVCF.tbi"
CALL_ARGS=()
LIST4="$TMPDIR/jasmine_4types.list"
: > "$LIST4"
for caller in delly manta smoove; do
  invcf=""
  case "$caller" in
    delly) invcf="$MED/results/delly_per_sample/$SID/$SID.delly.all.vcf.gz" ;;
    manta) invcf="$MED/results/manta_per_sample/$SID/results/variants/diploidSV.vcf.gz" ;;
    smoove) invcf="$MED/results/smoove_per_sample/$SID/$SID-smoove.genotyped.vcf.gz" ;;
  esac
  if [ -s "$invcf" ]; then
    raw="$NORMDIR/$SID.$caller.norm.raw.vcf"
    gz="$NORMDIR/$SID.$caller.norm.vcf.gz"
    jasmine_vcf="$TMPDIR/$SID.$caller.jasmine.4types.vcf"
    python3 "$MED/scripts/normalize_vcf.py" "$invcf" "$caller" "$SID" "$raw"
    "$BCFTOOLS" view -i 'INFO/SVTYPE!="BND" && INFO/SVTYPE!="TRA"' "$raw" -Oz -o "$gz" || true
    if [ -s "$gz" ]; then "$BCFTOOLS" index -t -f "$gz" || true; fi
    "$BCFTOOLS" view -i 'INFO/SVTYPE!="BND" && INFO/SVTYPE!="TRA"' "$raw" -Ov -o "$jasmine_vcf" || true
    if [ -s "$jasmine_vcf" ] && [ "$("$BCFTOOLS" view -H "$jasmine_vcf" | wc -l)" -gt 0 ]; then echo "$jasmine_vcf" >> "$LIST4"; fi
    bndraw="$NORMDIR/$SID.$caller.bnd.raw.vcf"
    bndgz="$NORMDIR/$SID.$caller.bnd.vcf.gz"
    "$BCFTOOLS" view -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="TRA"' "$raw" -Oz -o "$bndgz" || true
    if [ -s "$bndgz" ]; then "$BCFTOOLS" index -t -f "$bndgz" || true; CALL_ARGS+=("$caller=$bndgz"); fi
  fi
done
parts=()
if [ "$(wc -l < "$LIST4")" -ge 2 ]; then
  rawj="$TMPDIR/$SID.jasmine.4types.raw.vcf"
  fixedj="$TMPDIR/$SID.jasmine.4types.headerfix.vcf"
  sortj="$TMPDIR/$SID.jasmine.4types.sorted.vcf.gz"
  "$JAVA" -Xmx12g -jar "$JASMINE" file_list="$LIST4" out_file="$rawj" min_support=2 max_dist=500 threads=2 --ignore_strand --normalize_type > "$TMPDIR/jasmine.stdout" 2> "$TMPDIR/jasmine.stderr"
  if [ -s "$rawj" ] && [ "$(grep -vc '^#' "$rawj" || true)" -gt 0 ]; then
    awk '/^#CHROM/ {
      print "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Breakpoint strand orientation emitted by Jasmine\">"
      print "##INFO=<ID=CIGAR,Number=1,Type=String,Description=\"CIGAR string\">"
      print "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"Event identifier\">"
      print "##INFO=<ID=HOMSEQ,Number=1,Type=String,Description=\"Breakpoint homology sequence\">"
      print "##INFO=<ID=JUNCTION_QUAL,Number=1,Type=String,Description=\"Junction quality\">"
      print "##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description=\"Insertion length\">"
      print "##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description=\"Insertion sequence\">"
      print "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized phred-scaled genotype likelihoods\">"
      print "##FORMAT=<ID=PR,Number=2,Type=Integer,Description=\"Paired-read support\">"
      print "##FORMAT=<ID=SR,Number=2,Type=Integer,Description=\"Split-read support\">"
    } {print}' "$rawj" > "$fixedj"
    "$BCFTOOLS" sort -T "$TMPDIR/bcftools_sort_4" -Oz -o "$sortj" "$fixedj"
    "$BCFTOOLS" index -t -f "$sortj"
    parts+=("$sortj")
  fi
fi
bndraw="$TMPDIR/$SID.bnd.2caller.raw.vcf"
if [ "${#CALL_ARGS[@]}" -ge 2 ]; then
  python3 "$MED/scripts/merge_bnd_support.py" "$SID" "$bndraw" 1000 "${CALL_ARGS[@]}"
  if [ -s "$bndraw" ] && [ "$(grep -vc '^#' "$bndraw" || true)" -gt 0 ]; then
    bndgz="$TMPDIR/$SID.bnd.2caller.sorted.vcf.gz"
    "$BCFTOOLS" sort -T "$TMPDIR/bcftools_sort_bnd" -Oz -o "$bndgz" "$bndraw"
    "$BCFTOOLS" index -t -f "$bndgz"
    parts+=("$bndgz")
  fi
fi
if [ "${#parts[@]}" -eq 0 ]; then
  # create empty but valid VCF from reference header not trivial; write status and exit 0
  echo -e "sample_id\t$SID\nstatus\tno_2caller_consensus\ntime\t$(date '+%F %T %Z')" > "$DONE"
  echo "no consensus $SID"
  exit 0
fi
if [ "${#parts[@]}" -eq 1 ]; then
  cp -f "${parts[0]}" "$OUTVCF"
  cp -f "${parts[0]}.tbi" "$OUTVCF.tbi"
else
  "$BCFTOOLS" concat -a -Oz -o "$TMPDIR/$SID.concat.vcf.gz" "${parts[@]}"
  "$BCFTOOLS" sort -T "$TMPDIR/bcftools_sort_final" -Oz -o "$OUTVCF" "$TMPDIR/$SID.concat.vcf.gz"
  "$BCFTOOLS" index -t -f "$OUTVCF"
fi
{
  echo -e "sample_id\t$SID"
  echo -e "records\t$("$BCFTOOLS" view -H "$OUTVCF" | wc -l)"
  echo -e "delly_available\t$([ -s "$MED/results/delly_per_sample/$SID/$SID.delly.all.vcf.gz" ] && echo 1 || echo 0)"
  echo -e "manta_available\t$([ -s "$MED/results/manta_per_sample/$SID/results/variants/diploidSV.vcf.gz" ] && echo 1 || echo 0)"
  echo -e "smoove_available\t$([ -s "$MED/results/smoove_per_sample/$SID/$SID-smoove.genotyped.vcf.gz" ] && echo 1 || echo 0)"
  echo -e "time\t$(date '+%F %T %Z')"
} > "$DONE"
echo "consensus done $SID"
