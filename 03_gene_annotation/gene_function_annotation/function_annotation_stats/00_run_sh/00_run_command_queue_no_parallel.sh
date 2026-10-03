#!/usr/bin/env bash
set -euo pipefail

if [[ $# != 4 ]]; then
  echo "Usage: $0 JOBS COMMAND_FILE STDOUT_LOG STDERR_LOG" >&2
  exit 2
fi

jobs="$1"
command_file="$2"
stdout_log="$3"
stderr_log="$4"

mkdir -p "$(dirname "${stdout_log}")" "$(dirname "${stderr_log}")"

running=0
while IFS= read -r cmd || [[ -n "${cmd}" ]]; do
  [[ -z "${cmd}" ]] && continue
  {
    echo "[$(date)] START ${cmd}"
    bash -lc "${cmd}"
    echo "[$(date)] DONE ${cmd}"
  } >> "${stdout_log}" 2>> "${stderr_log}" &
  running=$((running + 1))
  if (( running >= jobs )); then
    wait -n
    running=$((running - 1))
  fi
done < "${command_file}"

wait
