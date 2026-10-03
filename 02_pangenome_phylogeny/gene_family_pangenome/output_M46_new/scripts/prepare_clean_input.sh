#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=path/to/project/1.orthology_family_2/data/pangenome_data/pep_18_genome_new
OUTPUT_ROOT=path/to/project/1.orthology_family_2/output_M46_new
CLEAN_DIR="$OUTPUT_ROOT/input_clean_18_genome"
QC_DIR="$OUTPUT_ROOT/qc"
MANIFEST="$QC_DIR/clean_input_manifest.tsv"

if [[ -e "$CLEAN_DIR" ]]; then
  printf 'Refusing to overwrite existing clean input directory: %s\n' "$CLEAN_DIR" >&2
  exit 1
fi

mkdir -p "$CLEAN_DIR" "$QC_DIR"
printf 'sample\toriginal_path\tclean_path\tsequences\tterminal_dots_removed\toriginal_sha256\tclean_sha256\n' > "$MANIFEST"

mapfile -t FASTA_FILES < <(find -L "$INPUT_DIR" -maxdepth 1 -type f -name '*.fa' | sort)
if [[ "${#FASTA_FILES[@]}" -ne 18 ]]; then
  printf 'Expected 18 readable FASTA files, found %s\n' "${#FASTA_FILES[@]}" >&2
  exit 1
fi

total_sequences=0
total_dots_removed=0

for original in "${FASTA_FILES[@]}"; do
  sample=$(basename "$original" .fa)
  resolved=$(readlink -f "$original")
  clean="$CLEAN_DIR/$sample.fa"

  awk '
    /^>/ { print; next }
    {
      gsub(/\./, "")
      print
    }
  ' "$resolved" > "$clean"

  sequences=$(grep -c '^>' "$clean")
  original_sequences=$(grep -c '^>' "$resolved")
  dots_removed=$(awk '!/^>/ { n += gsub(/\./, "") } END { print n+0 }' "$resolved")
  clean_dots=$(awk '!/^>/ { n += gsub(/\./, "") } END { print n+0 }' "$clean")
  empty_sequences=$(awk '
    /^>/ {
      if (seen_header && sequence_length == 0) empty++
      seen_header=1
      sequence_length=0
      next
    }
    {
      gsub(/[[:space:]]/, "")
      sequence_length += length($0)
    }
    END {
      if (seen_header && sequence_length == 0) empty++
      print empty+0
    }
  ' "$clean")

  if [[ "$sequences" -ne "$original_sequences" || "$clean_dots" -ne 0 || "$empty_sequences" -ne 0 ]]; then
    printf 'Clean-input validation failed for %s\n' "$sample" >&2
    exit 1
  fi

  original_sha256=$(sha256sum "$resolved" | cut -d ' ' -f 1)
  clean_sha256=$(sha256sum "$clean" | cut -d ' ' -f 1)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$resolved" "$clean" "$sequences" "$dots_removed" \
    "$original_sha256" "$clean_sha256" >> "$MANIFEST"

  total_sequences=$((total_sequences + sequences))
  total_dots_removed=$((total_dots_removed + dots_removed))
done

printf 'input_files=18\n' > "$QC_DIR/clean_input_summary.txt"
printf 'total_proteins=%s\n' "$total_sequences" >> "$QC_DIR/clean_input_summary.txt"
printf 'terminal_dots_removed=%s\n' "$total_dots_removed" >> "$QC_DIR/clean_input_summary.txt"
printf 'source_files_modified=0\n' >> "$QC_DIR/clean_input_summary.txt"
printf 'status=PASS\n' >> "$QC_DIR/clean_input_summary.txt"

touch "$CLEAN_DIR/.validated_clean_input"
cat "$QC_DIR/clean_input_summary.txt"
