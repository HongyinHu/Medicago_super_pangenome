#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 DONE_MARKER LOCK_FILE COMMAND..." >&2
  exit 2
fi

done_marker="$1"
lock_file="$2"
shift 2

if [[ -s "${done_marker}" ]]; then
  exit 0
fi

mkdir -p "$(dirname "${done_marker}")" "$(dirname "${lock_file}")"
exec 9>"${lock_file}"
if ! flock -n 9; then
  echo "[SKIP locked] ${lock_file}"
  exit 0
fi

if [[ -s "${done_marker}" ]]; then
  exit 0
fi

"$@"
date > "${done_marker}"
