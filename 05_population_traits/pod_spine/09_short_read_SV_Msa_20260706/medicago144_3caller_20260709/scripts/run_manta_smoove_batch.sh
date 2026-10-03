#!/usr/bin/env bash
set -euo pipefail
MED=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706/medicago144_3caller_20260709
BATCH_SIZE=${BATCH_SIZE:-4}
N=$(($(wc -l < "$MED/input/medicago144_sample_bams.tsv")-1))
NGROUPS=$(( (N + BATCH_SIZE - 1) / BATCH_SIZE ))
AID=${SLURM_ARRAY_TASK_ID:-${1:-}}
[ -n "$AID" ] || { echo "Need array id" >&2; exit 2; }
if [ "$AID" -le "$NGROUPS" ]; then
  export SLURM_ARRAY_TASK_ID="$AID"
  bash "$MED/scripts/run_batch.sh" manta
else
  export SLURM_ARRAY_TASK_ID=$((AID-NGROUPS))
  bash "$MED/scripts/run_batch.sh" smoove
fi
