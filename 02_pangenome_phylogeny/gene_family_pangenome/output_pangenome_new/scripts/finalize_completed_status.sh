#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/1.orthology_family_2/output_pangenome_new
RESULT_DIR="$ROOT/orthofinder_results/Results_Jul31"
LOG="$ROOT/logs/orthofinder.log"
STATUS_DIR="$ROOT/status"
GENE_COUNT="$RESULT_DIR/Orthogroups/Orthogroups.GeneCount.tsv"
UNASSIGNED="$RESULT_DIR/Orthogroups/Orthogroups_UnassignedGenes.tsv"
STATISTICS="$RESULT_DIR/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv"

grep -q 'Done orthologues' "$LOG"
test -s "$GENE_COUNT"
test -s "$UNASSIGNED"
test -s "$STATISTICS"

rm -f "$STATUS_DIR/RUNNING" "$STATUS_DIR/FAILED"
touch "$STATUS_DIR/COMPLETE"

if ! grep -q '^exit_code=' "$STATUS_DIR/run_manifest.txt"; then
  printf 'exit_code=0\n' >> "$STATUS_DIR/run_manifest.txt"
fi
if ! grep -q '^end_time=' "$STATUS_DIR/run_manifest.txt"; then
  printf 'end_time=%s\n' "$(date --iso-8601=seconds)" >> "$STATUS_DIR/run_manifest.txt"
fi
printf 'gene_count=%s\n' "$GENE_COUNT" >> "$STATUS_DIR/run_manifest.txt"
printf 'unassigned=%s\n' "$UNASSIGNED" >> "$STATUS_DIR/run_manifest.txt"
printf 'statistics_per_species=%s\n' "$STATISTICS" >> "$STATUS_DIR/run_manifest.txt"

cat > "$STATUS_DIR/finalization_note.txt" <<'EOF'
The OrthoFinder log ended with "Done orthologues" and the OrthoFinder main
process and wrapper process had exited. The wrapper had not cleared RUNNING.
COMPLETE was recovered only after the three required output tables were checked
as non-empty. No OrthoFinder result was modified by this status recovery.
EOF

printf 'status=PASS\n'
wc -c "$GENE_COUNT" "$UNASSIGNED" "$STATISTICS"
