#!/usr/bin/env bash
set -euo pipefail

STEP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN="path/to/project/N_3.call_SV/01.read_based_dualref_hifi"
SNIFFLES="$STEP/tools/sniffles2_skip_single"

task_file="${1:?task list required}"
task_id="${2:?task id required}"
CHUNK_THREADS="${CHUNK_THREADS:-1}"

line="$(sed -n "${task_id}p" "$task_file")"
if [[ -z "$line" ]]; then
  echo "No task line for task_id=$task_id in $task_file" >&2
  exit 2
fi

IFS=$'\t' read -r ref_alias ref_source ref_fa species chunk_id chunk_vcf chunk_bed <<< "$line"
bam="$RUN/03_per_sample/$ref_source/bam/$species.sorted.bam"
out_dir="$STEP/results/$ref_alias/genotype_chunks/$species"
mkdir -p "$out_dir" "$STEP/logs"
out_vcf="$out_dir/$chunk_id.vcf"
tmp_vcf="$out_dir/$chunk_id.tmp.vcf"
done_flag="$out_dir/$chunk_id.done"

if [[ -s "$out_vcf" && -s "$done_flag" ]]; then
  echo "SKIP done $ref_alias $species $chunk_id"
  exit 0
fi

rm -f "$tmp_vcf" "$done_flag"
"$SNIFFLES" \
  --input "$bam" \
  --genotype-vcf "$chunk_vcf" \
  --vcf "$tmp_vcf" \
  --reference "$ref_fa" \
  --regions "$chunk_bed" \
  --threads "$CHUNK_THREADS" \
  > "$STEP/logs/genotype_chunk.$ref_alias.$species.$chunk_id.out" \
  2> "$STEP/logs/genotype_chunk.$ref_alias.$species.$chunk_id.err"

test -s "$tmp_vcf"
mv "$tmp_vcf" "$out_vcf"
touch "$done_flag"
echo "DONE $ref_alias $species $chunk_id"
