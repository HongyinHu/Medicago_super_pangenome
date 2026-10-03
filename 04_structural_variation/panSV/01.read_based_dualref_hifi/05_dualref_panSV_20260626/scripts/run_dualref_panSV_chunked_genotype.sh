#!/usr/bin/env bash
set -euo pipefail

STEP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN="path/to/project/N_3.call_SV/01.read_based_dualref_hifi"
export STEP RUN TMPDIR="$STEP/tmp"
mkdir -p "$STEP"/{logs,tmp,summary,config/lists,results}

PY=path/to/home/anaconda3/envs/biosofeware/bin/python
BGZIP=path/to/home/anaconda3/envs/biosofeware/bin/bgzip
TABIX=path/to/home/anaconda3/envs/biosofeware/bin/tabix
SNIFFLES="$STEP/tools/sniffles2_skip_single"
MINIMAP2=path/to/home/anaconda3/envs/minimap2_env/bin/minimap2
ALIGN_THREADS=${ALIGN_THREADS:-24}
CHUNK_MAX_RECORDS=${CHUNK_MAX_RECORDS:-5000}
CHUNK_PADDING=${CHUNK_PADDING:-20000}
CHUNK_JOBS=${CHUNK_JOBS:-4}
CHUNK_THREADS=${CHUNK_THREADS:-1}

log() { echo "[$(date '+%F %T')] $*" | tee -a "$STEP/logs/run_dualref_panSV_chunked_genotype.log"; }

need_tool() {
  local x="$1"
  [[ -x "$x" ]] || { echo "ERROR: tool not executable: $x" >&2; exit 2; }
}

for x in "$PY" "$BGZIP" "$TABIX" "$SNIFFLES" "$MINIMAP2"; do
  need_tool "$x"
done

compress_index() {
  local vcf="$1"
  local tmp="$vcf.gz.inprogress"
  rm -f "$tmp" "$vcf.gz" "$vcf.gz.tbi"
  "$BGZIP" -f -c "$vcf" > "$tmp" || { rm -f "$tmp"; return 1; }
  mv "$tmp" "$vcf.gz"
  "$TABIX" -f -p vcf "$vcf.gz"
}
export -f compress_index
export BGZIP TABIX

genotype_chunk_one() {
  local ref_alias="$1" ref_source="$2" ref_fa="$3" species="$4" chunk_id="$5" chunk_vcf="$6" chunk_bed="$7"
  local bam out_dir tmp_vcf out_vcf done_flag
  bam="$RUN/03_per_sample/$ref_source/bam/$species.sorted.bam"
  out_dir="$STEP/results/$ref_alias/genotype_chunks/$species"
  mkdir -p "$out_dir"
  out_vcf="$out_dir/$chunk_id.vcf"
  tmp_vcf="$out_dir/$chunk_id.tmp.vcf"
  done_flag="$out_dir/$chunk_id.done"
  if [[ -s "$out_vcf" && -s "$done_flag" ]]; then
    return 0
  fi
  rm -f "$tmp_vcf" "$done_flag"
  "$SNIFFLES" \
    --input "$bam" \
    --genotype-vcf "$chunk_vcf" \
    --vcf "$tmp_vcf" \
    --reference "$ref_fa" \
    --regions "$chunk_bed" \
    --threads "$CHUNK_THREADS" \
    > "$STEP/logs/genotype_chunk.$ref_alias.$species.$chunk_id.out" \
    2> "$STEP/logs/genotype_chunk.$ref_alias.$species.$chunk_id.err"
  mv "$tmp_vcf" "$out_vcf"
  touch "$done_flag"
}
export -f genotype_chunk_one
export STEP RUN SNIFFLES CHUNK_THREADS

log "checking required discovery outputs"
while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
  [[ "$ref_alias" == "ref_alias" ]] && continue
  [[ -s "$STEP/results/$ref_alias/panSV/panSV.$ref_alias.discovery.vcf" ]] || { echo "missing discovery VCF for $ref_alias" >&2; exit 2; }
  [[ -s "$ref_fa.fai" ]] || { echo "missing fai: $ref_fa.fai" >&2; exit 2; }
done < "$STEP/config/refs.tsv"

log "splitting discovery VCFs into genotype chunks"
while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
  [[ "$ref_alias" == "ref_alias" ]] && continue
  chunk_dir="$STEP/results/$ref_alias/genotype_chunks/chunks"
  manifest="$STEP/results/$ref_alias/genotype_chunks/chunks.tsv"
  if [[ ! -s "$manifest" ]]; then
    rm -rf "$chunk_dir"
    mkdir -p "$chunk_dir"
    "$PY" "$STEP/scripts/split_discovery_vcf.py" \
      --vcf "$STEP/results/$ref_alias/panSV/panSV.$ref_alias.discovery.vcf" \
      --fai "$ref_fa.fai" \
      --ref-alias "$ref_alias" \
      --out-dir "$chunk_dir" \
      --manifest "$manifest" \
      --max-records "$CHUNK_MAX_RECORDS" \
      --padding "$CHUNK_PADDING" \
      > "$STEP/logs/split_discovery.$ref_alias.out" \
      2> "$STEP/logs/split_discovery.$ref_alias.err"
  fi
done < "$STEP/config/refs.tsv"
touch "$STEP/summary/genotype_chunks_split.done"

log "running chunked Sniffles2 genotype with RefA/RefB interleaved by chunk"
chunk_args="$STEP/tmp/genotype_chunks.args"
: > "$chunk_args"
while read -r species; do
  max_chunks="$(awk -F'\t' 'FNR>1 && (FNR-1)>m {m=FNR-1} END {print m+0}' "$STEP"/results/*/genotype_chunks/chunks.tsv)"
  for idx in $(seq 1 "$max_chunks"); do
    while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
      [[ "$ref_alias" == "ref_alias" ]] && continue
      manifest="$STEP/results/$ref_alias/genotype_chunks/chunks.tsv"
      row="$(awk -F'\t' -v n="$idx" 'NR==n+1 {print; exit}' "$manifest")"
      [[ -n "$row" ]] || continue
      IFS=$'\t' read -r chunk_id chrom start end count chunk_vcf chunk_bed <<< "$row"
      out_dir="$STEP/results/$ref_alias/genotype_chunks/$species"
      if [[ -s "$out_dir/$chunk_id.vcf" && -s "$out_dir/$chunk_id.done" ]]; then
        continue
      fi
      echo -e "$ref_alias\t$ref_source\t$ref_fa\t$species\t$chunk_id\t$chunk_vcf\t$chunk_bed" >> "$chunk_args"
    done < "$STEP/config/refs.tsv"
  done
done < "$STEP/config/species.list"
if [[ -s "$chunk_args" ]]; then
  xargs -P "$CHUNK_JOBS" -n 7 bash -c 'set -euo pipefail; genotype_chunk_one "$@"' _ < "$chunk_args"
fi
touch "$STEP/summary/chunked_genotype.done"

log "concatenating chunked genotype VCFs"
while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
  [[ "$ref_alias" == "ref_alias" ]] && continue
  manifest="$STEP/results/$ref_alias/genotype_chunks/chunks.tsv"
  while read -r species; do
    raw="$STEP/results/$ref_alias/genotype/$species.$ref_alias.panSV.gt.raw.vcf"
    sorted="$STEP/results/$ref_alias/genotype/$species.$ref_alias.panSV.gt.vcf"
    if [[ ! -s "$sorted.gz" || ! -s "$sorted.gz.tbi" ]]; then
      "$PY" "$STEP/scripts/concat_chunk_vcfs.py" \
        --manifest "$manifest" \
        --chunk-dir "$STEP/results/$ref_alias/genotype_chunks/$species" \
        --species "$species" \
        --ref-alias "$ref_alias" \
        --out-vcf "$raw" \
        > "$STEP/logs/concat_chunks.$ref_alias.$species.out" \
        2> "$STEP/logs/concat_chunks.$ref_alias.$species.err"
      "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$raw" --out-vcf "$sorted" --fai "$ref_fa.fai"
      compress_index "$sorted"
    fi
  done < "$STEP/config/species.list"
done < "$STEP/config/refs.tsv"
touch "$STEP/summary/genotype.done"

log "merging per-species genotype VCFs and making PAV matrices"
while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
  [[ "$ref_alias" == "ref_alias" ]] && continue
  list="$STEP/config/lists/$ref_alias/18species.$ref_alias.genotype.list"
  : > "$list"
  while read -r species; do
    echo -e "$species\t$STEP/results/$ref_alias/genotype/$species.$ref_alias.panSV.gt.vcf.gz" >> "$list"
  done < "$STEP/config/species.list"
  raw="$STEP/results/$ref_alias/panSV/panSV.$ref_alias.genotyped.18species.raw.vcf"
  sorted="$STEP/results/$ref_alias/panSV/panSV.$ref_alias.genotyped.18species.vcf"
  "$PY" "$STEP/scripts/merge_genotyped_vcfs.py" --list "$list" --out "$raw" \
    > "$STEP/logs/merge_genotyped.$ref_alias.out" \
    2> "$STEP/logs/merge_genotyped.$ref_alias.err"
  "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$raw" --out-vcf "$sorted" --fai "$ref_fa.fai"
  compress_index "$sorted"
  "$PY" "$STEP/scripts/make_pav_matrix.py" \
    --vcf "$sorted.gz" \
    --matrix "$STEP/results/$ref_alias/panSV/panSV.$ref_alias.PAV.matrix.tsv" \
    --stats "$STEP/results/$ref_alias/panSV/panSV.$ref_alias.stats.tsv"
done < "$STEP/config/refs.tsv"
touch "$STEP/summary/reference_specific_panSV.done"

log "building preliminary cross-reference meta-panSV table"
mkdir -p "$STEP/results/meta-panSV"
refa_fa="$(awk -F'\t' '$1=="RefA"{print $3}' "$STEP/config/refs.tsv")"
refb_fa="$(awk -F'\t' '$1=="RefB"{print $3}' "$STEP/config/refs.tsv")"
paf="$STEP/results/meta-panSV/RefA_to_RefB.minimap2.paf"
if [[ ! -s "$paf" ]]; then
  "$MINIMAP2" -x asm5 -t "$ALIGN_THREADS" "$refb_fa" "$refa_fa" > "$paf" 2> "$STEP/logs/minimap2.RefA_to_RefB.err"
fi
"$PY" "$STEP/scripts/cross_ref_match.py" \
  --refa-vcf "$STEP/results/RefA/panSV/panSV.RefA.genotyped.18species.vcf.gz" \
  --refb-vcf "$STEP/results/RefB/panSV/panSV.RefB.genotyped.18species.vcf.gz" \
  --paf "$paf" \
  --out-events "$STEP/results/meta-panSV/meta-panSV.events.tsv" \
  --out-evidence "$STEP/results/meta-panSV/meta-panSV.matching.evidence.tsv" \
  > "$STEP/logs/cross_ref_match.out" \
  2> "$STEP/logs/cross_ref_match.err"
touch "$STEP/summary/meta_panSV.done"

date > "$STEP/summary/dualref_panSV.done"
log "done"
