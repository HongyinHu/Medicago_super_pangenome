#!/usr/bin/env bash
set -euo pipefail

RUN="$1"
MAX_JOBS="$2"
THREADS="$3"
READS_PREFIX_ALLOW="${4:-}"

TASKS="${RUN}/01_metadata/readbased_tasks.tsv"
LOCK_DIR="${RUN}/locks"
STATUS_DIR="${RUN}/status"
SCRIPT="${RUN}/scripts/run_readbased_task.sh"
HOST=$(hostname)
WORKER_LOG="${RUN}/logs/workers/${HOST}.worker.$(date +%Y%m%d_%H%M%S).log"

mkdir -p "${LOCK_DIR}" "${STATUS_DIR}" "${RUN}/logs/workers"
exec > >(tee -a "${WORKER_LOG}") 2>&1

echo "WORKER_START host=${HOST} max_jobs=${MAX_JOBS} threads=${THREADS} reads_prefix_allow=${READS_PREFIX_ALLOW:-ALL} time=$(date -Is)"

count_active() {
  jobs -pr | wc -l
}

all_claimed_or_done() {
  local pending=0
  while IFS=$'\t' read -r task_id rest; do
    [[ "${task_id}" == "task_id" ]] && continue
    if [[ ! -e "${STATUS_DIR}/${task_id}.done" && ! -e "${STATUS_DIR}/${task_id}.failed" && ! -d "${LOCK_DIR}/${task_id}.lock" ]]; then
      pending=$((pending + 1))
      break
    fi
  done < "${TASKS}"
  [[ "${pending}" -eq 0 ]]
}

while true; do
  made_claim=0
  while (( $(count_active) < MAX_JOBS )); do
    claimed=0
    while IFS= read -r line; do
      task_id=$(printf "%s" "${line}" | cut -f1)
      [[ "${task_id}" == "task_id" ]] && continue
      if [[ -n "${READS_PREFIX_ALLOW}" ]]; then
        reads_path=$(printf "%s" "${line}" | cut -f7)
        [[ "${reads_path}" == "${READS_PREFIX_ALLOW}"* ]] || continue
      fi
      [[ -e "${STATUS_DIR}/${task_id}.done" ]] && continue
      [[ -e "${STATUS_DIR}/${task_id}.failed" ]] && continue
      lock="${LOCK_DIR}/${task_id}.lock"
      if mkdir "${lock}" 2>/dev/null; then
        {
          echo "host=${HOST}"
          echo "pid=$$"
          echo "threads=${THREADS}"
          echo "start=$(date -Is)"
          echo "line=${line}"
        } > "${lock}/claim.txt"
        echo "CLAIM task=${task_id} host=${HOST} time=$(date -Is)"
        bash "${SCRIPT}" "${RUN}" "${line}" "${THREADS}" &
        claimed=1
        made_claim=1
        break
      fi
    done < "${TASKS}"
    if [[ "${claimed}" -eq 0 ]]; then
      break
    fi
  done

  if (( $(count_active) == 0 )); then
    if all_claimed_or_done; then
      echo "WORKER_NO_PENDING host=${HOST} time=$(date -Is)"
      break
    fi
  fi

  sleep 30
done

wait || true
echo "WORKER_END host=${HOST} time=$(date -Is)"
