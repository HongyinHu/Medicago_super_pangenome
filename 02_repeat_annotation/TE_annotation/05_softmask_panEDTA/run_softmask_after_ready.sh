#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/N_1.EDTA_single
PAN="$BASE/03_panEDTA"
LAI="$BASE/04_LAI_panEDTA"
SOFT="$BASE/05_softmask_panEDTA"
LOGDIR="$SOFT/logs"
export PATH=path/to/home/anaconda3/envs/EDTA_env/bin:$PATH
GENOMES=(
  genome_395.fa genome_410.fa genome_436.fa genome_454.fa genome_457.fa
  genome_461.fa genome_468.fa genome_472.fa genome_474a.fa genome_474b.fa
  genome_482.fa genome_A17.fa genome_M22.fa genome_M46.fa genome_Mar.fa
  genome_Mpo.fa genome_Mru.fa genome_Msa1.fa genome_Msa2.fa genome_R108.fa
  genome_ZM4.fa
)

mkdir -p "$SOFT" "$LOGDIR" "$SOFT/bed" "$SOFT/fasta"
exec >>"$LOGDIR/finalizer.log" 2>&1

ts() {
  date '+%F %T %Z'
}

structural_done_count() {
  (find "$PAN/logs/structural_resume/reconciled" -maxdepth 1 -type f -name '*.done' 2>/dev/null || true) | wc -l
}

structural_failed_count() {
  (find "$PAN/logs/structural_resume/failed" -maxdepth 1 -type f -name '*.failed' 2>/dev/null || true) | wc -l
}

lai_ready() {
  [ -s "$LAI/LAI_summary.tsv" ] && return 0
  local done na failed
  done=$( (find "$LAI/logs/recover_lai/done" -maxdepth 1 -type f -name '*.done' 2>/dev/null || true) | wc -l )
  na=$( (find "$LAI/logs/recover_lai/not_applicable" -maxdepth 1 -type f -name '*.done' 2>/dev/null || true) | wc -l )
  failed=$( (find "$LAI/logs/recover_lai/failed" -maxdepth 1 -type f -name '*.failed' 2>/dev/null || true) | wc -l )
  [ "$failed" -eq 0 ] && [ $((done + na)) -eq 21 ]
}

pick_existing_masked() {
  local g=$1
  local mod="${g}.mod"
  local c
  for c in \
    "$PAN/${mod}.panEDTA.softmasked.fa" \
    "$PAN/${mod}.panEDTA.masked" \
    "$PAN/${mod}.masked" \
    "$PAN/${mod}.EDTA.anno/${mod}.masked" \
    "$PAN/${mod}.EDTA.raw/${mod}.masked"
  do
    [ -s "$c" ] && {
      printf '%s\n' "$c"
      return 0
    }
  done
  return 1
}

out_to_bed() {
  local out=$1
  local bed=$2
  awk 'BEGIN{OFS="\t"} $1 ~ /^[0-9]+$/ && $6 ~ /^[0-9]+$/ && $7 ~ /^[0-9]+$/ {
    start=$6-1;
    if (start < 0) start=0;
    print $5, start, $7
  }' "$out" | sort -k1,1 -k2,2n | bedtools merge -i - > "$bed"
}

make_softmask_one() {
  local g=$1
  local mod="${g}.mod"
  local sample="${g%.fa}"
  local genome="$PAN/$mod"
  local rmout="$PAN/${mod}.panEDTA.out"
  local bed="$SOFT/bed/${sample}.panEDTA.repeatmasker.bed"
  local fa="$SOFT/fasta/${sample}.panEDTA.softmasked.fa"
  local src

  if [ -s "$fa" ]; then
    echo "SKIP $(ts) $g existing=$fa"
    return 0
  fi

  if src=$(pick_existing_masked "$g"); then
    cp -p "$src" "$fa"
    echo "COPY $(ts) $g src=$src out=$fa"
    return 0
  fi

  [ -s "$genome" ] || {
    echo "ERROR $(ts) $g missing_genome=$genome"
    return 1
  }
  [ -s "$rmout" ] || {
    echo "ERROR $(ts) $g missing_rmout=$rmout"
    return 1
  }

  out_to_bed "$rmout" "$bed.tmp"
  mv "$bed.tmp" "$bed"
  bedtools maskfasta -soft -fi "$genome" -bed "$bed" -fo "$fa.tmp"
  mv "$fa.tmp" "$fa"
  echo "MASK $(ts) $g bed=$bed out=$fa"
}

write_manifest() {
  local manifest="$SOFT/manifest.tsv"
  printf 'sample\tgenome_mod\trmout\tsoftmask_fasta\tbed\tsoftmask_size_bytes\n' > "$manifest.tmp"
  for g in "${GENOMES[@]}"; do
    local sample="${g%.fa}"
    local fa="$SOFT/fasta/${sample}.panEDTA.softmasked.fa"
    local bed="$SOFT/bed/${sample}.panEDTA.repeatmasker.bed"
    local size=0
    [ -s "$fa" ] && size=$(stat -c '%s' "$fa")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$sample" "$PAN/${g}.mod" "$PAN/${g}.mod.panEDTA.out" "$fa" "$bed" "$size" >> "$manifest.tmp"
  done
  mv "$manifest.tmp" "$manifest"
  find "$SOFT/fasta" -maxdepth 1 -type f -name '*.panEDTA.softmasked.fa' | wc -l > "$SOFT/softmask.count"
}

run_softmask_all() {
  command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR $(ts) bedtools_not_found"
    return 1
  }

  echo "SOFTMASK_START $(ts)"
  local failed=0
  for g in "${GENOMES[@]}"; do
    make_softmask_one "$g" || failed=$((failed + 1))
  done
  write_manifest
  local count
  count=$(cat "$SOFT/softmask.count")
  echo "SOFTMASK_DONE $(ts) count=$count failed=$failed"
  [ "$count" -eq 21 ] && [ "$failed" -eq 0 ]
}

echo "FINALIZER_START $(ts) host=$(hostname)"
while true; do
  sdone=$(structural_done_count)
  sfail=$(structural_failed_count)
  soft_count=$(find "$SOFT/fasta" -maxdepth 1 -type f -name '*.panEDTA.softmasked.fa' 2>/dev/null | wc -l || true)
  echo "CHECK $(ts) structural_done=$sdone/21 structural_failed=$sfail softmask_count=$soft_count"
  if [ "$sfail" -gt 0 ]; then
    echo "STOP $(ts) structural_failed=$sfail"
    exit 2
  fi
  if [ "$sdone" -eq 21 ]; then
    run_softmask_all
    if lai_ready; then
      echo "ALL_READY $(ts) softmask=21 lai_ready=yes"
    else
      echo "SOFTMASK_READY_WAIT_LAI $(ts)"
    fi
    exit 0
  fi
  sleep 120
done
