#!/usr/bin/env bash
set -euo pipefail
MED=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706/medicago144_3caller_20260709
CALLER=${1:?caller required}
BATCH_SIZE=${BATCH_SIZE:-4}
N=$(($(wc -l < "$MED/input/medicago144_sample_bams.tsv")-1))
GROUP=${SLURM_ARRAY_TASK_ID:-${2:-}}
[ -n "$GROUP" ] || { echo "Need array group id" >&2; exit 2; }
start=$(( (GROUP - 1) * BATCH_SIZE + 1 ))
end=$(( GROUP * BATCH_SIZE ))
[ "$end" -gt "$N" ] && end="$N"
echo "[$(date '+%F %T')] batch caller=$CALLER group=$GROUP samples=${start}-${end} host=$(hostname)"
for task in $(seq "$start" "$end"); do
  case "$CALLER" in
    delly) env -u SLURM_ARRAY_TASK_ID bash "$MED/scripts/run_delly_one.sh" "$task" ;;
    manta) env -u SLURM_ARRAY_TASK_ID bash "$MED/scripts/run_manta_one.sh" "$task" ;;
    smoove) env -u SLURM_ARRAY_TASK_ID bash "$MED/scripts/run_smoove_one.sh" "$task" ;;
    consensus) env -u SLURM_ARRAY_TASK_ID bash "$MED/scripts/run_consensus_one.sh" "$task" ;;
    *) echo "unknown caller $CALLER" >&2; exit 2 ;;
  esac
done
echo "[$(date '+%F %T')] batch done caller=$CALLER group=$GROUP"
