#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/N_1.coding_gene_anno
PI="$BASE/04_transcript_prediction/PASA_inputs"
PASA="$BASE/04_transcript_prediction/PASA"
LOCK="$PI/locks/collect_all_pasa_transcripts.lock"
PROV=${PROV:?PROV is required}
NODE=$(hostname)

mkdir -p "$PI/locks" "$PROV"
if find "$PASA/locks" -mindepth 1 -maxdepth 1 -type d -name '*.lock' -print -quit 2>/dev/null | grep -q .; then
  echo "[ERROR] PASA sample lock exists" >&2
  exit 2
fi
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "[ERROR] collector lock exists: $LOCK" >&2
  exit 3
fi
cleanup() {
  rm -f -- "$LOCK/owner"
  rmdir "$LOCK" 2>/dev/null || true
}
trap cleanup EXIT
printf 'host=%s\npid=%s\nstarted=%s\n' "$NODE" "$$" "$(date '+%F %T %Z')" > "$LOCK/owner"

echo "[START] $(date '+%F %T %Z') node=$NODE"
bash "$PI/run/collect_all_pasa_transcripts.sh"

expected=$(awk 'NR > 1 {n++} END {print n+0}' "$BASE/01_genome_versions/genome_manifest.tsv")
actual=$(awk 'NR > 1 {n++} END {print n+0}' "$PI/manifests/pasa_transcripts_manifest.tsv")
not_ok=$(awk -F'\t' 'NR > 1 && $7 != "ok" {n++} END {print n+0}' "$PI/manifests/pasa_transcripts_manifest.tsv")
duplicates=$(awk -F'\t' 'NR > 1 {n[$1]++} END {for (s in n) if (n[s] > 1) d++; print d+0}' "$PI/manifests/pasa_transcripts_manifest.tsv")

msa1=$(awk -F'\t' '$1 == "genome_Msa1" {print; exit}' "$PI/manifests/pasa_transcripts_manifest.tsv")
mru=$(awk -F'\t' '$1 == "genome_Mru" {print; exit}' "$PI/manifests/pasa_transcripts_manifest.tsv")
msa1_pb=$(printf '%s\n' "$msa1" | awk -F'\t' '{print $5}')
msa1_existing=$(printf '%s\n' "$msa1" | awk -F'\t' '{print $6}')
mru_total=$(printf '%s\n' "$mru" | awk -F'\t' '{print $3}')
mru_existing=$(printf '%s\n' "$mru" | awk -F'\t' '{print $6}')

repeated_prefix=$(grep '^>' "$PI/combined/genome_Msa1.pasa_transcripts.fa" \
  | awk -F'|' '$1 == ">genome_Msa1" && $4 == "genome_Msa1" && $2 == $5 {n++} END {print n+0}')

[ "$actual" -eq "$expected" ]
[ "$not_ok" -eq 0 ]
[ "$duplicates" -eq 0 ]
[ "$msa1_pb" -eq 59429 ]
[ "$msa1_existing" -eq 0 ]
[ "$mru_total" -eq 38448 ]
[ "$mru_existing" -eq 38458 ]
[ "$repeated_prefix" -eq 0 ]

cp -p "$PI/manifests/pasa_transcripts_manifest.tsv" "$PROV/new_pasa_transcripts_manifest.tsv"
sha256sum "$PI/combined/genome_Msa1.pasa_transcripts.fa" > "$PROV/new_combined_sha256.txt"
{
  echo "finished=$(date '+%F %T %Z')"
  echo "expected_samples=$expected"
  echo "actual_samples=$actual"
  echo "not_ok=$not_ok"
  echo "duplicates=$duplicates"
  echo "msa1_pacbio=$msa1_pb"
  echo "msa1_existing=$msa1_existing"
  echo "mru_total=$mru_total"
  echo "mru_existing=$mru_existing"
  echo "repeated_prefix=$repeated_prefix"
} > "$PROV/collector_validation.txt"
echo "[DONE] $(date '+%F %T %Z') samples=$actual msa1_pacbio=$msa1_pb"
