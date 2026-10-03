#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_5.pod_coil/03_SNP_INDEL_SV_GWAS_144_coil_ge1_20260714
SOURCE=path/to/project/N_4.pod_spiny/28_GWAS_spiny/05_sv_genotyping
REF=path/to/project/38.medicago_resequence/2.call_SNP_new/00_reference_index/ref_genome.T2T.fa
BAM=path/to/project/38.medicago_resequence/2.call_SNP_new/04_markdup_all/bam/qh3.dedup.bam
BAI=${BAM%.bam}.bai
SITES="$SOURCE/candidate_sites/medicago143.sites.filtered.nsamp2.smoove.vcf.gz"
SITES_INDEXED="$SOURCE/candidate_sites/medicago143.sites.filtered.nsamp2.vcf.gz"
SMOOVE=path/to/home/anaconda3/envs/smoove/bin/smoove
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
export PATH=path/to/home/anaconda3/envs/smoove/bin:path/to/home/anaconda3/envs/panpop/bin:/usr/bin:/bin

OUTDIR="$ROOT/03_sv_gwas/qh3_genotype"
OUT="$OUTDIR/qh3-smoove.genotyped.vcf.gz"
DONE="$OUTDIR/qh3.done"
mkdir -p "$OUTDIR" "$ROOT/logs" "$ROOT/summary"
LOG="$ROOT/logs/run_qh3_sv_genotype.log"
exec > >(tee -a "$LOG") 2>&1

test -s "$BAM"
test -s "$BAI"
test -s "$SITES"
EXPECTED=$("$BCFTOOLS" index -n "$SITES_INDEXED")
if [[ ! -s "$OUT" || ! -s "$OUT.csi" ]] \
    || ! "$BCFTOOLS" view -h "$OUT" >/dev/null 2>&1 \
    || [[ "$("$BCFTOOLS" index -n "$OUT" 2>/dev/null || true)" != "$EXPECTED" ]] \
    || [[ "$("$BCFTOOLS" query -l "$OUT" 2>/dev/null || true)" != qh3 ]]; then
  rm -f "$DONE"
  "$SMOOVE" genotype --name qh3 --outdir "$OUTDIR" --fasta "$REF" --removepr \
    --processes "${SLURM_CPUS_PER_TASK:-4}" --vcf "$SITES" "$BAM"
fi
test -s "$OUT"
test -s "$OUT.csi"
[[ "$("$BCFTOOLS" index -n "$OUT")" == "$EXPECTED" ]]
[[ "$("$BCFTOOLS" query -l "$OUT")" == qh3 ]]
printf 'completed\t%s\nhost\t%s\nsample\tqh3\nrecords\t%s\n' \
  "$(date '+%F %T %Z')" "$(hostname)" "$EXPECTED" > "$DONE"
cp "$DONE" "$ROOT/summary/qh3_sv_genotype.done"
