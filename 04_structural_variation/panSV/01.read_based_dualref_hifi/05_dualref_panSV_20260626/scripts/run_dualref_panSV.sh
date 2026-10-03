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
JAVA=path/to/home/anaconda3/envs/panpop/bin/java
JASMINE_JAR="$STEP/tools/jasmine.jar"

MIN_SVLEN=50
SVTYPES=DEL,INS,DUP,INV
FILTER_JOBS=${FILTER_JOBS:-12}
CONSENSUS_JOBS=${CONSENSUS_JOBS:-8}
GENOTYPE_JOBS=${GENOTYPE_JOBS:-1}
GENOTYPE_THREADS=${GENOTYPE_THREADS:-4}
ALIGN_THREADS=${ALIGN_THREADS:-32}

log() { echo "[$(date '+%F %T')] $*" | tee -a "$STEP/logs/run_dualref_panSV.log"; }

need_tool() {
  local x="$1"
  [[ -x "$x" ]] || { echo "ERROR: tool not executable: $x" >&2; exit 2; }
}

for x in "$PY" "$BGZIP" "$TABIX" "$SNIFFLES" "$MINIMAP2" "$JAVA"; do
  need_tool "$x"
done
[[ -s "$JASMINE_JAR" ]] || { echo "ERROR: Jasmine jar not found: $JASMINE_JAR" >&2; exit 2; }

read_refs() {
  tail -n +2 "$STEP/config/refs.tsv"
}

caller_suffix() {
  case "$1" in
    pbsv) echo "pbsv pbsv" ;;
    sniffles2) echo "sniffles2 sniffles2" ;;
    cutesv) echo "cutesv cutesv" ;;
    *) echo "ERROR unknown caller $1" >&2; return 2 ;;
  esac
}

compress_index() {
  local vcf="$1"
  local tmp="$vcf.gz.inprogress"
  rm -f "$tmp" "$vcf.gz" "$vcf.gz.tbi"
  "$BGZIP" -f -c "$vcf" > "$tmp" || { rm -f "$tmp"; return 1; }
  mv "$tmp" "$vcf.gz"
  "$TABIX" -f -p vcf "$vcf.gz"
}

filter_one() {
  local ref_alias="$1" ref_source="$2" ref_fa="$3" species="$4" caller="$5"
  local pair subdir suffix in_vcf out sorted
  pair="$(caller_suffix "$caller")"
  subdir="${pair%% *}"
  suffix="${pair##* }"
  in_vcf="$RUN/03_per_sample/$ref_source/$subdir/$species.$suffix.vcf.gz"
  out="$STEP/results/$ref_alias/filtered/$species.$ref_alias.$caller.filtered.vcf"
  sorted="$STEP/results/$ref_alias/filtered/$species.$ref_alias.$caller.filtered.sorted.vcf"
  if [[ -s "$sorted.gz" && -s "$sorted.gz.tbi" ]]; then
    return 0
  fi
  mkdir -p "$(dirname "$out")"
  "$PY" "$STEP/scripts/vcf_filter.py" \
    --in-vcf "$in_vcf" --out-vcf "$out" --sample "$species" --caller "$caller" \
    --ref-fai "$ref_fa.fai" --min-svlen "$MIN_SVLEN" --svtypes "$SVTYPES" \
    > "$STEP/logs/filter.$ref_alias.$species.$caller.out" \
    2> "$STEP/logs/filter.$ref_alias.$species.$caller.err"
  "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$out" --out-vcf "$sorted" --fai "$ref_fa.fai"
  compress_index "$sorted"
}
export -f filter_one caller_suffix compress_index
export PY BGZIP TABIX STEP RUN MIN_SVLEN SVTYPES

consensus_one() {
  local ref_alias="$1" ref_source="$2" ref_fa="$3" species="$4"
  local dir list raw sorted
  dir="$STEP/results/$ref_alias/consensus"
  mkdir -p "$dir" "$STEP/config/lists/$ref_alias"
  list="$STEP/config/lists/$ref_alias/$species.$ref_alias.callers.list"
  raw="$dir/$species.$ref_alias.consensus.jasmine.raw.vcf"
  sorted="$dir/$species.$ref_alias.consensus.vcf"
  if [[ -s "$sorted.gz" && -s "$sorted.gz.tbi" ]]; then
    return 0
  fi
  : > "$list"
  for caller in pbsv sniffles2 cutesv; do
    echo "$STEP/results/$ref_alias/filtered/$species.$ref_alias.$caller.filtered.sorted.vcf" >> "$list"
  done
  "$JAVA" -jar "$JASMINE_JAR" \
    file_list="$list" \
    out_file="$raw" \
    min_support=2 \
    max_dist=500 \
    threads=4 \
    --nonlinear_dist \
    --ignore_strand \
    --normalize_type \
    > "$STEP/logs/consensus.$ref_alias.$species.out" \
    2> "$STEP/logs/consensus.$ref_alias.$species.err"
  "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$raw" --out-vcf "$sorted" --fai "$ref_fa.fai"
  compress_index "$sorted"
}
export -f consensus_one
export JAVA JASMINE_JAR

genotype_one() {
  local ref_alias="$1" ref_source="$2" ref_fa="$3" species="$4"
  local bam disc dir raw sorted
  bam="$RUN/03_per_sample/$ref_source/bam/$species.sorted.bam"
  disc="$STEP/results/$ref_alias/panSV/panSV.$ref_alias.discovery.vcf"
  dir="$STEP/results/$ref_alias/genotype"
  mkdir -p "$dir"
  raw="$dir/$species.$ref_alias.panSV.gt.raw.vcf"
  sorted="$dir/$species.$ref_alias.panSV.gt.vcf"
  if [[ -s "$sorted.gz" && -s "$sorted.gz.tbi" ]]; then
    return 0
  fi
  "$SNIFFLES" --input "$bam" --genotype-vcf "$disc" --vcf "$raw" \
    --reference "$ref_fa" --threads "$GENOTYPE_THREADS" \
    > "$STEP/logs/genotype.$ref_alias.$species.out" \
    2> "$STEP/logs/genotype.$ref_alias.$species.err"
  "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$raw" --out-vcf "$sorted" --fai "$ref_fa.fai"
  compress_index "$sorted"
}
export -f genotype_one
export SNIFFLES GENOTYPE_THREADS

log "checking inputs"
"$PY" "$STEP/scripts/check_inputs.py" --step "$STEP" --run-dir "$RUN" --out "$STEP/config/input_manifest.tsv"

log "filtering caller VCFs"
filter_cmds="$STEP/tmp/filter.args"
: > "$filter_cmds"
while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
  [[ "$ref_alias" == "ref_alias" ]] && continue
  while read -r species; do
    for caller in pbsv sniffles2 cutesv; do
      echo "$ref_alias $ref_source $ref_fa $species $caller" >> "$filter_cmds"
    done
  done < "$STEP/config/species.list"
done < "$STEP/config/refs.tsv"
xargs -P "$FILTER_JOBS" -n 5 bash -c 'set -euo pipefail; filter_one "$@"' _ < "$filter_cmds"
touch "$STEP/summary/filter.done"

log "building within-species consensus VCFs"
cons_cmds="$STEP/tmp/consensus.args"
: > "$cons_cmds"
while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
  [[ "$ref_alias" == "ref_alias" ]] && continue
  while read -r species; do
    echo "$ref_alias $ref_source $ref_fa $species" >> "$cons_cmds"
  done < "$STEP/config/species.list"
done < "$STEP/config/refs.tsv"
xargs -P "$CONSENSUS_JOBS" -n 4 bash -c 'set -euo pipefail; consensus_one "$@"' _ < "$cons_cmds"
touch "$STEP/summary/consensus.done"

log "building reference-specific pan-SV discovery VCFs"
while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
  [[ "$ref_alias" == "ref_alias" ]] && continue
  mkdir -p "$STEP/results/$ref_alias/panSV" "$STEP/config/lists/$ref_alias"
  list="$STEP/config/lists/$ref_alias/18species.$ref_alias.consensus.list"
  : > "$list"
  while read -r species; do
    echo "$STEP/results/$ref_alias/consensus/$species.$ref_alias.consensus.vcf" >> "$list"
  done < "$STEP/config/species.list"
  raw="$STEP/results/$ref_alias/panSV/panSV.$ref_alias.discovery.jasmine.raw.vcf"
  sorted="$STEP/results/$ref_alias/panSV/panSV.$ref_alias.discovery.vcf"
  if [[ ! -s "$sorted.gz" || ! -s "$sorted.gz.tbi" ]]; then
    "$JAVA" -jar "$JASMINE_JAR" \
      file_list="$list" \
      out_file="$raw" \
      min_support=1 \
      max_dist=1000 \
      threads=32 \
      --nonlinear_dist \
      --ignore_strand \
      --normalize_type \
      > "$STEP/logs/discovery.$ref_alias.out" \
      2> "$STEP/logs/discovery.$ref_alias.err"
    "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$raw" --out-vcf "$sorted" --fai "$ref_fa.fai"
    compress_index "$sorted"
  fi
done < "$STEP/config/refs.tsv"
touch "$STEP/summary/discovery.done"

log "genotyping reference-specific pan-SV candidates with Sniffles2"
gt_cmds="$STEP/tmp/genotype.args"
: > "$gt_cmds"
while read -r species; do
  while IFS=$'\t' read -r ref_alias ref_source ref_fa; do
    [[ "$ref_alias" == "ref_alias" ]] && continue
    echo "$ref_alias $ref_source $ref_fa $species" >> "$gt_cmds"
  done < "$STEP/config/refs.tsv"
done < "$STEP/config/species.list"
xargs -P "$GENOTYPE_JOBS" -n 4 bash -c 'set -euo pipefail; genotype_one "$@"' _ < "$gt_cmds"
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
