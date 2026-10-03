#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=path/to/project/1.orthology_family_2/data/pangenome_data/pep_18_genome
OUTPUT_ROOT=path/to/project/1.orthology_family_2/output_pangenome_new
REPORT="$OUTPUT_ROOT/qc/sequence_character_audit.tsv"

mkdir -p "$OUTPUT_ROOT/qc"
printf 'sample\tsequences\tdot_characters\tsequences_ending_dot\tsequences_with_internal_dot\tstar_characters\tdash_characters\tnon_letter_other\n' > "$REPORT"

for fasta in "$INPUT_DIR"/*.fa; do
  sample=$(basename "$fasta" .fa)
  awk -v sample="$sample" '
    function finalize_sequence(    i, c, internal_dot) {
      if (!have_sequence) return
      sequences++
      internal_dot=0
      for (i=1; i<=length(sequence); i++) {
        c=substr(sequence, i, 1)
        if (c == ".") {
          dots++
          if (i < length(sequence)) internal_dot=1
        } else if (c == "*") {
          stars++
        } else if (c == "-") {
          dashes++
        } else if (c !~ /[A-Za-z]/) {
          other++
        }
      }
      if (substr(sequence, length(sequence), 1) == ".") terminal_dot_sequences++
      if (internal_dot) internal_dot_sequences++
    }
    /^>/ {
      finalize_sequence()
      sequence=""
      have_sequence=1
      next
    }
    {
      gsub(/[[:space:]]/, "")
      sequence=sequence $0
    }
    END {
      finalize_sequence()
      printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", sample, sequences, dots, terminal_dot_sequences, internal_dot_sequences, stars, dashes, other
    }
  ' "$fasta" >> "$REPORT"
done

cat "$REPORT"
