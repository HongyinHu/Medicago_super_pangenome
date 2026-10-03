#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=path/to/project/1.orthology_family_2/data/pangenome_data/pep_18_genome_new
OUTPUT_ROOT=path/to/project/1.orthology_family_2/output_M46_new
QC_DIR="$OUTPUT_ROOT/qc"

mkdir -p "$QC_DIR"
MANIFEST="$QC_DIR/input_manifest.tsv"

mapfile -t FASTA_FILES < <(find -L "$INPUT_DIR" -maxdepth 1 -type f -name '*.fa' | sort)
if [[ "${#FASTA_FILES[@]}" -ne 18 ]]; then
  printf 'Expected 18 readable FASTA files, found %s\n' "${#FASTA_FILES[@]}" >&2
  exit 1
fi

printf 'sample\tsymlink\tresolved_path\tsize_bytes\tmtime\tsequences\tduplicate_primary_ids\tempty_sequences\n' > "$MANIFEST"

total_sequences=0
for fasta in "${FASTA_FILES[@]}"; do
  sample=$(basename "$fasta" .fa)
  resolved=$(readlink -f "$fasta")
  if [[ ! -s "$resolved" ]]; then
    printf 'Missing or empty FASTA target: %s -> %s\n' "$fasta" "$resolved" >&2
    exit 1
  fi

  size_bytes=$(stat -Lc '%s' "$resolved")
  mtime=$(stat -Lc '%y' "$resolved")
  sequences=$(grep -c '^>' "$resolved")
  duplicate_ids=$(awk '
    /^>/ {
      id=$1
      sub(/^>/, "", id)
      seen[id]++
      if (seen[id] == 2) duplicates++
    }
    END { print duplicates+0 }
  ' "$resolved")
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
  ' "$resolved")

  if [[ "$sequences" -eq 0 || "$duplicate_ids" -ne 0 || "$empty_sequences" -ne 0 ]]; then
    printf 'FASTA QC failed for %s: sequences=%s duplicates=%s empty=%s\n' \
      "$sample" "$sequences" "$duplicate_ids" "$empty_sequences" >&2
    exit 1
  fi

  total_sequences=$((total_sequences + sequences))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$fasta" "$resolved" "$size_bytes" "$mtime" \
    "$sequences" "$duplicate_ids" "$empty_sequences" >> "$MANIFEST"
done

printf 'input_files=18\n' > "$QC_DIR/input_qc_summary.txt"
printf 'total_proteins=%s\n' "$total_sequences" >> "$QC_DIR/input_qc_summary.txt"
printf 'duplicate_primary_ids=0\n' >> "$QC_DIR/input_qc_summary.txt"
printf 'empty_sequences=0\n' >> "$QC_DIR/input_qc_summary.txt"
printf 'status=PASS\n' >> "$QC_DIR/input_qc_summary.txt"

cat "$QC_DIR/input_qc_summary.txt"
