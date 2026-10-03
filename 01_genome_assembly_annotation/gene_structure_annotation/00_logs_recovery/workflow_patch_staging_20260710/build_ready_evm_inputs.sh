#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/N_1.coding_gene_anno
E="$BASE/05_EVM_integration"
PASA_MANIFEST="$BASE/04_transcript_prediction/PASA_inputs/manifests/pasa_transcripts_manifest.tsv"
GEMOMA_ROOT="$BASE/03_homology_prediction/GeMoMa_strict_v2"
mkdir -p "$E/inputs" "$E/status" "$E/logs" "$E/locks"
READY="$E/ready_samples.tsv"
: > "$READY"
printf 'sample\tgenome\tbraker_gff3\tbraker_source\tgemoma_gff\tpasa_gff3\tstatus\tinput_signature\n' > "$READY"

input_signature() {
  local path
  for path in "$@"; do
    [ -s "$path" ] || return 1
    stat -Lc '%n\t%s\t%Y' "$path"
  done | sha256sum | awk '{print $1}'
}

for genome in "$BASE"/01_genome_versions/unmask/genome_*.unmasked.fa; do
  [ -s "$genome" ] || continue
  sample=$(basename "$genome" .unmasked.fa)
  pasa_done="$BASE/04_transcript_prediction/PASA/$sample/pasa.done"
  pasa_gff="$BASE/04_transcript_prediction/PASA/$sample/${sample}.sqlite.pasa_assemblies.gff3"
  gemoma_done="$GEMOMA_ROOT/$sample/gemoma_strict_v2.done"
  gemoma_gff="$GEMOMA_ROOT/$sample/final_annotation.gff"
  illumina_count=$(awk -F'\t' -v s="$sample" 'NR>1 && $1==s{print $4+0; found=1; exit} END{if(!found) print 0}' "$PASA_MANIFEST")
  braker_gff=""; braker_src=""
  if [ "$illumina_count" -gt 0 ]; then
    braker_gff="$BASE/02_abinitio_prediction/BRAKER3_RNA_protein/$sample/braker.gff3"; braker_src="BRAKER3_RNA_protein"
    braker_done="$BASE/02_abinitio_prediction/BRAKER3_RNA_protein/$sample/braker3_rna_protein.done"
  else
    braker_gff="$BASE/02_abinitio_prediction/BRAKER3_protein/$sample/braker.gff3"; braker_src="BRAKER3_protein"
    braker_done="$BASE/02_abinitio_prediction/BRAKER3_protein/$sample/braker3_protein.done"
  fi
  missing=()
  [ -s "$pasa_done" ] || missing+=(pasa_done)
  [ -s "$pasa_gff" ] || missing+=(pasa_gff3)
  [ -s "$gemoma_done" ] || missing+=(gemoma_done)
  [ -s "$gemoma_gff" ] || missing+=(gemoma_gff)
  [ -s "$braker_done" ] || missing+=(braker_done)
  [ -s "$braker_gff" ] || missing+=(braker_gff3)
  if [ ${#missing[@]} -eq 0 ]; then
    status=ready
    signature=$(input_signature "$genome" "$braker_gff" "$gemoma_gff" "$pasa_gff")
    out="$E/inputs/$sample"
    mkdir -p "$out"
    ln -sfn "$genome" "$out/genome.fa"
    ln -sfn "$braker_gff" "$out/braker.gff3.raw"
    ln -sfn "$gemoma_gff" "$out/gemoma.gff3.raw"
    ln -sfn "$pasa_gff" "$out/pasa_assemblies.gff3.raw"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$genome" "$braker_gff" "$braker_src" "$gemoma_gff" "$pasa_gff" "$status" "$signature" >> "$READY"
  else
    status="missing:${missing[*]}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\n' "$sample" "$genome" "${braker_gff:-NA}" "${braker_src:-NA}" "$gemoma_gff" "$pasa_gff" "$status" >> "$READY"
  fi
done
awk -F'\t' 'NR>1{n[$7]++} END{for(k in n) print k,n[k]}' "$READY" | sort > "$E/ready_summary.txt"
