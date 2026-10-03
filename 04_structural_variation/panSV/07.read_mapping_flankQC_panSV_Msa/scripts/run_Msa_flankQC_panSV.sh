#!/usr/bin/env bash
set -euo pipefail

BASE="path/to/project/N_3.call_SV"
RUN="$BASE/01.read_based_dualref_hifi"
PREV="$RUN/05_dualref_panSV_20260626"
STEP="$BASE/07.read_mapping_flankQC_panSV_Msa"
REF_ALIAS="Msa"
REF_SOURCE="Msa"
REF_FA="$RUN/02_ref_prepare/Msa/ref.fa"
SVG_DIR="$BASE/02.assembly_svgap_dualref/Msa_ref/CombinedSV"

export STEP RUN REF_FA TMPDIR="$STEP/tmp"
mkdir -p "$STEP"/{scripts,tools,logs,tmp,summary,config/lists,results/Msa/{filtered,consensus,consensus_flankQC,panSV,qc}}

PY=${PY:-path/to/home/anaconda3/envs/biosofeware/bin/python}
BGZIP=${BGZIP:-path/to/home/anaconda3/envs/biosofeware/bin/bgzip}
TABIX=${TABIX:-path/to/home/anaconda3/envs/biosofeware/bin/tabix}
SAMTOOLS=${SAMTOOLS:-path/to/home/anaconda3/envs/biosofeware/bin/samtools}
MOSDEPTH=${MOSDEPTH:-path/to/home/anaconda3/envs/minimap2_env/bin/mosdepth}
JAVA=${JAVA:-path/to/home/anaconda3/envs/panpop/bin/java}
JASMINE_JAR="$STEP/tools/jasmine.jar"

MIN_SVLEN=${MIN_SVLEN:-50}
SVTYPES=${SVTYPES:-DEL,INS,DUP,INV,TRA}
FILTER_JOBS=${FILTER_JOBS:-12}
CONSENSUS_JOBS=${CONSENSUS_JOBS:-6}
FLANK_QC_JOBS=${FLANK_QC_JOBS:-6}
CONSENSUS_THREADS=${CONSENSUS_THREADS:-4}
DISCOVERY_THREADS=${DISCOVERY_THREADS:-32}
FLANK_SIZE=${FLANK_SIZE:-500}
MOSDEPTH_THREADS=${MOSDEPTH_THREADS:-2}
MIN_EACH_FLANK_DEPTH=${MIN_EACH_FLANK_DEPTH:-5}
MIN_MEAN_FLANK_DEPTH=${MIN_MEAN_FLANK_DEPTH:-8}

log() { echo "[$(date '+%F %T')] $*" | tee -a "$STEP/logs/run_Msa_flankQC_panSV.log"; }

need_file() {
  [[ -s "$1" ]] || { echo "ERROR: missing file: $1" >&2; exit 2; }
}

need_exec() {
  [[ -x "$1" ]] || { echo "ERROR: missing executable: $1" >&2; exit 2; }
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
  local species="$1" caller="$2"
  local pair subdir suffix in_vcf out sorted
  pair="$(caller_suffix "$caller")"
  subdir="${pair%% *}"
  suffix="${pair##* }"
  in_vcf="$RUN/03_per_sample/$REF_SOURCE/$subdir/$species.$suffix.vcf.gz"
  out="$STEP/results/Msa/filtered/$species.Msa.$caller.filtered.vcf"
  sorted="$STEP/results/Msa/filtered/$species.Msa.$caller.filtered.sorted.vcf"
  if [[ -s "$sorted.gz" && -s "$sorted.gz.tbi" ]]; then
    return 0
  fi
  need_file "$in_vcf"
  "$PY" "$STEP/scripts/vcf_filter_msa_with_tra.py" \
    --in-vcf "$in_vcf" --out-vcf "$out" --sample "$species" --caller "$caller" \
    --ref-fai "$REF_FA.fai" --min-svlen "$MIN_SVLEN" --svtypes "$SVTYPES" \
    --interchrom-tra-only \
    > "$STEP/logs/filter.Msa.$species.$caller.out" \
    2> "$STEP/logs/filter.Msa.$species.$caller.err"
  "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$out" --out-vcf "$sorted" --fai "$REF_FA.fai"
  compress_index "$sorted"
}
export -f filter_one caller_suffix compress_index need_file
export PY BGZIP TABIX STEP RUN REF_SOURCE REF_FA MIN_SVLEN SVTYPES

consensus_one() {
  local species="$1"
  local dir list raw sorted
  dir="$STEP/results/Msa/consensus"
  mkdir -p "$dir" "$STEP/config/lists/Msa"
  list="$STEP/config/lists/Msa/$species.Msa.callers.list"
  raw="$dir/$species.Msa.consensus.jasmine.raw.vcf"
  sorted="$dir/$species.Msa.consensus.vcf"
  if [[ -s "$sorted.gz" && -s "$sorted.gz.tbi" ]]; then
    return 0
  fi
  : > "$list"
  for caller in pbsv sniffles2 cutesv; do
    echo "$STEP/results/Msa/filtered/$species.Msa.$caller.filtered.sorted.vcf" >> "$list"
  done
  "$JAVA" -jar "$JASMINE_JAR" \
    file_list="$list" \
    out_file="$raw" \
    min_support=2 \
    max_dist=500 \
    threads="$CONSENSUS_THREADS" \
    --nonlinear_dist \
    --ignore_strand \
    --normalize_type \
    > "$STEP/logs/consensus.Msa.$species.out" \
    2> "$STEP/logs/consensus.Msa.$species.err"
  "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$raw" --out-vcf "$sorted" --fai "$REF_FA.fai"
  compress_index "$sorted"
}
export -f consensus_one
export JAVA JASMINE_JAR CONSENSUS_THREADS

flank_qc_one() {
  local species="$1"
  local in_vcf bam out qc tmp_prefix sorted
  in_vcf="$STEP/results/Msa/consensus/$species.Msa.consensus.vcf"
  bam="$RUN/03_per_sample/$REF_SOURCE/bam/$species.sorted.bam"
  out="$STEP/results/Msa/consensus_flankQC/$species.Msa.consensus.flankQC.vcf"
  sorted="$STEP/results/Msa/consensus_flankQC/$species.Msa.consensus.flankQC.sorted.vcf"
  qc="$STEP/results/Msa/qc/$species.Msa.consensus.flankQC.tsv"
  tmp_prefix="$STEP/tmp/flankQC.$species"
  if [[ -s "$sorted.gz" && -s "$sorted.gz.tbi" && -s "$qc" ]]; then
    return 0
  fi
  need_file "$in_vcf"
  need_file "$bam"
  "$PY" "$STEP/scripts/flank_coverage_filter_vcf.py" \
    --in-vcf "$in_vcf" \
    --out-vcf "$out" \
    --qc-tsv "$qc" \
    --bam "$bam" \
    --ref-fai "$REF_FA.fai" \
    --samtools "$SAMTOOLS" \
    --mosdepth "$MOSDEPTH" \
    --threads "$MOSDEPTH_THREADS" \
    --tmp-prefix "$tmp_prefix" \
    --flank "$FLANK_SIZE" \
    --min-each-flank-depth "$MIN_EACH_FLANK_DEPTH" \
    --min-mean-flank-depth "$MIN_MEAN_FLANK_DEPTH" \
    > "$STEP/logs/flankQC.Msa.$species.out" \
    2> "$STEP/logs/flankQC.Msa.$species.err"
  "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$out" --out-vcf "$sorted" --fai "$REF_FA.fai"
  compress_index "$sorted"
}
export -f flank_qc_one
export SAMTOOLS MOSDEPTH MOSDEPTH_THREADS FLANK_SIZE MIN_EACH_FLANK_DEPTH MIN_MEAN_FLANK_DEPTH

main() {
  log "checking tools and inputs"
  need_exec "$PY"
  need_exec "$BGZIP"
  need_exec "$TABIX"
  need_exec "$SAMTOOLS"
  need_exec "$MOSDEPTH"
  need_exec "$JAVA"
  need_file "$REF_FA"
  need_file "$REF_FA.fai"
  need_file "$JASMINE_JAR"
  need_file "$STEP/config/species.list"
  need_file "$SVG_DIR/All.DELs.50bplarge.bed.combined.sorted.txt"

  log "filtering per-caller VCFs with SVTYPE=$SVTYPES"
  filter_cmds="$STEP/tmp/filter.args"
  : > "$filter_cmds"
  while read -r species; do
    [[ -z "$species" ]] && continue
    for caller in pbsv sniffles2 cutesv; do
      echo "$species $caller" >> "$filter_cmds"
    done
  done < "$STEP/config/species.list"
  xargs -P "$FILTER_JOBS" -n 2 bash -c 'set -euo pipefail; filter_one "$@"' _ < "$filter_cmds"
  touch "$STEP/summary/filter.done"

  log "building within-species Jasmine consensus"
  cons_cmds="$STEP/tmp/consensus.args"
  : > "$cons_cmds"
  while read -r species; do
    [[ -z "$species" ]] && continue
    echo "$species" >> "$cons_cmds"
  done < "$STEP/config/species.list"
  xargs -P "$CONSENSUS_JOBS" -n 1 bash -c 'set -euo pipefail; consensus_one "$@"' _ < "$cons_cmds"
  touch "$STEP/summary/consensus.done"

  log "filtering species consensus VCFs by flank read coverage"
  qc_cmds="$STEP/tmp/flankQC.args"
  : > "$qc_cmds"
  while read -r species; do
    [[ -z "$species" ]] && continue
    echo "$species" >> "$qc_cmds"
  done < "$STEP/config/species.list"
  xargs -P "$FLANK_QC_JOBS" -n 1 bash -c 'set -euo pipefail; flank_qc_one "$@"' _ < "$qc_cmds"
  touch "$STEP/summary/flankQC.done"

  log "building Msa panSV discovery from flankQC consensus"
  mkdir -p "$STEP/results/Msa/panSV" "$STEP/config/lists/Msa"
  list="$STEP/config/lists/Msa/18species.Msa.flankQC.consensus.list"
  : > "$list"
  while read -r species; do
    [[ -z "$species" ]] && continue
    echo "$STEP/results/Msa/consensus_flankQC/$species.Msa.consensus.flankQC.sorted.vcf" >> "$list"
  done < "$STEP/config/species.list"
  raw="$STEP/results/Msa/panSV/panSV.Msa.flankQC.discovery.jasmine.raw.vcf"
  sorted="$STEP/results/Msa/panSV/panSV.Msa.flankQC.discovery.vcf"
  if [[ ! -s "$sorted.gz" || ! -s "$sorted.gz.tbi" ]]; then
    "$JAVA" -jar "$JASMINE_JAR" \
      file_list="$list" \
      out_file="$raw" \
      min_support=1 \
      max_dist=1000 \
      threads="$DISCOVERY_THREADS" \
      --nonlinear_dist \
      --ignore_strand \
      --normalize_type \
      > "$STEP/logs/discovery.Msa.out" \
      2> "$STEP/logs/discovery.Msa.err"
    "$PY" "$STEP/scripts/sort_vcf.py" --in-vcf "$raw" --out-vcf "$sorted" --fai "$REF_FA.fai"
    compress_index "$sorted"
  fi
  touch "$STEP/summary/discovery.done"

  log "making read-mapping PAV matrix from Jasmine support vectors"
  "$PY" "$STEP/scripts/make_pav_from_jasmine_vcf.py" \
    --vcf "$sorted.gz" \
    --species-list "$STEP/config/species.list" \
    --matrix "$STEP/results/Msa/panSV/panSV.Msa.flankQC.read_mapping.PAV.matrix.tsv" \
    --stats "$STEP/results/Msa/panSV/panSV.Msa.flankQC.read_mapping.stats.tsv" \
    > "$STEP/logs/make_pav.Msa.out" \
    2> "$STEP/logs/make_pav.Msa.err"
  touch "$STEP/summary/read_mapping_pav.done"

  log "integrating SVGAP as read-primary evidence only"
  "$PY" "$STEP/scripts/integrate_svgap_read_primary_msa.py" \
    --read-pav "$STEP/results/Msa/panSV/panSV.Msa.flankQC.read_mapping.PAV.matrix.tsv" \
    --svgap-dir "$SVG_DIR" \
    --out-events "$STEP/results/Msa/panSV/panSV.Msa.flankQC.read_primary.events.tsv" \
    --out-pav "$STEP/results/Msa/panSV/panSV.Msa.flankQC.read_primary.PAV.matrix.tsv" \
    --out-evidence "$STEP/results/Msa/panSV/panSV.Msa.flankQC.svgap_overlap_evidence.tsv" \
    --out-stats "$STEP/results/Msa/panSV/panSV.Msa.flankQC.read_primary.stats.tsv" \
    > "$STEP/logs/integrate_svgap.Msa.out" \
    2> "$STEP/logs/integrate_svgap.Msa.err"
  touch "$STEP/summary/svgap_evidence_integrated.done"

  date > "$STEP/summary/Msa_flankQC_panSV.done"
  log "done"
}

main "$@"
