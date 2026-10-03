#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=path/to/project/1.orthology_family_2/output_M46_new/input_clean_18_genome
OUTPUT_ROOT=path/to/project/1.orthology_family_2/output_M46_new
RESULT_DIR="$OUTPUT_ROOT/orthofinder_results"
LOG_DIR="$OUTPUT_ROOT/logs"
STATUS_DIR="$OUTPUT_ROOT/status"
CONDA_SH=path/to/home/anaconda3/etc/profile.d/conda.sh
THREADS=30

mkdir -p "$LOG_DIR" "$STATUS_DIR"
printf '%s\n' "$$" > "$STATUS_DIR/pid"

if [[ -e "$RESULT_DIR" ]]; then
  printf 'Refusing to overwrite existing result directory: %s\n' "$RESULT_DIR" >&2
  exit 1
fi

set +u
source "$CONDA_SH"
conda activate biosofeware
set -u

ORTHOFINDER=$(command -v orthofinder)
DIAMOND=$(command -v diamond)
MCL=$(command -v mcl)

{
  printf 'start_time=%s\n' "$(date --iso-8601=seconds)"
  printf 'host=%s\n' "$(hostname)"
  printf 'input_dir=%s\n' "$INPUT_DIR"
  printf 'result_dir=%s\n' "$RESULT_DIR"
  printf 'orthofinder=%s\n' "$ORTHOFINDER"
  printf 'diamond=%s\n' "$DIAMOND"
  printf 'mcl=%s\n' "$MCL"
  printf 'command=orthofinder -f %s -o %s -M msa -os -S diamond -t %s -a %s\n' "$INPUT_DIR" "$RESULT_DIR" "$THREADS" "$THREADS"
} > "$STATUS_DIR/run_manifest.txt"

"$ORTHOFINDER" -h > "$STATUS_DIR/orthofinder_help.txt" 2>&1
head -n 2 "$STATUS_DIR/orthofinder_help.txt" >> "$STATUS_DIR/run_manifest.txt"
"$DIAMOND" version >> "$STATUS_DIR/run_manifest.txt" 2>&1
"$MCL" --version >> "$STATUS_DIR/run_manifest.txt" 2>&1 || true

touch "$STATUS_DIR/RUNNING"
rm -f "$STATUS_DIR/FAILED" "$STATUS_DIR/COMPLETE"

set +e
"$ORTHOFINDER" \
  -f "$INPUT_DIR" \
  -o "$RESULT_DIR" \
  -M msa \
  -os \
  -S diamond \
  -t "$THREADS" \
  -a "$THREADS" \
  > "$LOG_DIR/orthofinder.log" 2>&1
exit_code=$?
set -e

rm -f "$STATUS_DIR/RUNNING"
printf 'exit_code=%s\n' "$exit_code" >> "$STATUS_DIR/run_manifest.txt"
printf 'end_time=%s\n' "$(date --iso-8601=seconds)" >> "$STATUS_DIR/run_manifest.txt"

if [[ "$exit_code" -ne 0 ]]; then
  touch "$STATUS_DIR/FAILED"
  exit "$exit_code"
fi

GENE_COUNT=$(find "$RESULT_DIR" -type f -path '*/Orthogroups/Orthogroups.GeneCount.tsv' | head -n 1)
UNASSIGNED=$(find "$RESULT_DIR" -type f -path '*/Orthogroups/Orthogroups_UnassignedGenes.tsv' | head -n 1)
STATS=$(find "$RESULT_DIR" -type f -path '*/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv' | head -n 1)

if [[ -z "$GENE_COUNT" || ! -s "$GENE_COUNT" || -z "$UNASSIGNED" || ! -s "$UNASSIGNED" || -z "$STATS" || ! -s "$STATS" ]]; then
  printf 'Required OrthoFinder output is missing after exit code 0\n' >&2
  touch "$STATUS_DIR/FAILED"
  exit 2
fi

printf 'gene_count=%s\n' "$GENE_COUNT" >> "$STATUS_DIR/run_manifest.txt"
printf 'unassigned=%s\n' "$UNASSIGNED" >> "$STATUS_DIR/run_manifest.txt"
printf 'statistics_per_species=%s\n' "$STATS" >> "$STATUS_DIR/run_manifest.txt"
touch "$STATUS_DIR/COMPLETE"
