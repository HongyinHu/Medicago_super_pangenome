#!/usr/bin/env bash
set -euo pipefail
sample=${1:?sample}
threads=${2:-8}
BASE=path/to/project/N_1.coding_gene_anno
E="$BASE/05_EVM_integration"
mkdir -p "$E/status" "$E/locks" "$E/logs" "$E/work"
lock="$E/locks/$sample.lock"
if mkdir "$lock" 2>/dev/null; then
  echo "host=$(hostname)" > "$lock/owner"
  echo "pid=$$" >> "$lock/owner"
  echo "started=$(date '+%F %T %Z')" >> "$lock/owner"
else
  echo "locked $sample" >&2; exit 9
fi
cleanup() {
  rc=$?
  if [ "$rc" -eq 0 ]; then
    rm -rf "$lock"
  else
    echo "failed rc=$rc $(date '+%F %T %Z') host=$(hostname)" > "$E/status/$sample.failed"
    echo "failed_rc=$rc $(date '+%F %T %Z')" >> "$lock/owner" 2>/dev/null || true
  fi
  exit "$rc"
}
trap cleanup EXIT
rm -f "$E/status/$sample.failed"
source ~/.bashrc 2>/dev/null || true
conda activate evm_env
export PATH="$CONDA_PREFIX/bin/EvmUtils:$PATH"
"$E/run/prepare_evm_sample.sh" "$sample"
out="$E/work/$sample"
data="$out/00_inputs"
run="$out/01_evm_run"
final="$out/02_final"
mkdir -p "$run" "$final"
cd "$run"
ln -sfn "$data/genome.fa" genome.fa
ln -sfn "$data/gene_predictions.gff3" gene_predictions.gff3
ln -sfn "$data/transcript_alignments.gff3" transcript_alignments.gff3
ln -sfn "$data/weights.txt" weights.txt
rm -rf partitions partitions_list.out commands.list evm_command_status.tsv evm_commands.stdout.log evm_commands.stderr.log
partition_EVM_inputs.pl --partition_dir partitions --genome genome.fa --gene_predictions gene_predictions.gff3 --transcript_alignments transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
write_EVM_commands.pl --genome genome.fa --weights "$PWD/weights.txt" --gene_predictions gene_predictions.gff3 --transcript_alignments transcript_alignments.gff3 --output_file_name evm.out --partitions partitions_list.out > commands.list
: > evm_command_status.tsv
: > evm_commands.stdout.log
: > evm_commands.stderr.log
fail=0
running=0
while IFS= read -r cmd; do
  [ -n "$cmd" ] || continue
  (
    bash -lc "$cmd"
    rc=$?
    printf '%s\t%s\n' "$rc" "$cmd" >> evm_command_status.tsv
    exit "$rc"
  ) >> evm_commands.stdout.log 2>> evm_commands.stderr.log &
  running=$((running+1))
  if [ "$running" -ge "$threads" ]; then
    wait -n || fail=1
    running=$((running-1))
  fi
done < commands.list
while [ "$running" -gt 0 ]; do
  wait -n || fail=1
  running=$((running-1))
done
if [ "$fail" -ne 0 ]; then
  echo "one_or_more_EVM_partition_commands_failed" >&2
  exit 20
fi
recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output_file_name evm.out --genome genome.fa
find partitions -name evm.out.gff3 -type f -print0 | sort -z | xargs -0 cat > "$final/$sample.evm.gff3"
/usr/bin/perl path/to/home/sofeware/PASApipeline-v2.3.3/misc_utilities/gff3_file_to_proteins.pl "$final/$sample.evm.gff3" genome.fa prot > "$final/$sample.evm.pep.fa"
/usr/bin/perl path/to/home/sofeware/PASApipeline-v2.3.3/misc_utilities/gff3_file_to_proteins.pl "$final/$sample.evm.gff3" genome.fa CDS > "$final/$sample.evm.cds.fa"
printf 'sample\tgenes\tpep\tcds\n%s\t%s\t%s\t%s\n' "$sample" "$(awk '$3=="gene"{n++} END{print n+0}' "$final/$sample.evm.gff3")" "$(grep -c '^>' "$final/$sample.evm.pep.fa" || true)" "$(grep -c '^>' "$final/$sample.evm.cds.fa" || true)" > "$final/$sample.evm.stats.tsv"
[ -s "$final/$sample.evm.gff3" ] && [ -s "$final/$sample.evm.pep.fa" ] && [ -s "$final/$sample.evm.cds.fa" ]
genes=$(awk '$3=="gene"{n++} END{print n+0}' "$final/$sample.evm.gff3")
pep=$(grep -c '^>' "$final/$sample.evm.pep.fa" || true)
cds=$(grep -c '^>' "$final/$sample.evm.cds.fa" || true)
signature=$(awk -F'\t' '$1=="input_signature"{print $2}' "$data/provenance.tsv")
[ "$genes" -gt 0 ] && [ "$pep" -gt 0 ] && [ "$cds" -gt 0 ] || { echo "empty_final_counts genes=$genes pep=$pep cds=$cds" >&2; exit 21; }
bad_aa_lines=$(grep -v '^>' "$final/$sample.evm.pep.fa" | LC_ALL=C grep -c '[^A-Z*]' || true)
[ "$bad_aa_lines" -eq 0 ] || { echo "invalid_protein_characters lines=$bad_aa_lines" >&2; exit 22; }
echo "done $(date '+%F %T %Z') signature=$signature genes=$genes pep=$pep cds=$cds" > "$E/status/$sample.done"
