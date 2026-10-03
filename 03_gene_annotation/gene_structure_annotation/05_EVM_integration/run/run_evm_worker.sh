#!/usr/bin/env bash
set -euo pipefail
node_label=${1:?node_label}
threads=${2:-8}
shift 2
BASE=path/to/project/N_1.coding_gene_anno
E="$BASE/05_EVM_integration"
mkdir -p "$E/logs" "$E/status"
echo "[$(date '+%F %T %Z')] START worker node=$node_label host=$(hostname) threads=$threads samples=$*"
for sample in "$@"; do
  if [ -s "$E/status/$sample.done" ]; then
    echo "[$(date '+%F %T %Z')] SKIP done $sample"
    continue
  fi
  if [ -d "$E/locks/$sample.lock" ]; then
    echo "[$(date '+%F %T %Z')] SKIP locked $sample"
    continue
  fi
  echo "[$(date '+%F %T %Z')] RUN $sample"
  if bash "$E/run/run_evm_sample.sh" "$sample" "$threads" > "$E/logs/${sample}.evm.stdout.log" 2> "$E/logs/${sample}.evm.stderr.log"; then
    echo "[$(date '+%F %T %Z')] DONE $sample"
  else
    rc=$?
    echo "[$(date '+%F %T %Z')] FAILED $sample rc=$rc"
    echo "failed rc=$rc $(date '+%F %T %Z')" > "$E/status/$sample.failed"
  fi
done
echo "[$(date '+%F %T %Z')] FINISH worker node=$node_label host=$(hostname)"
