#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project/N_4.pod_spiny/28_GWAS_spiny
STATUS="$ROOT/03_population_structure/summary/driver.status.tsv"
mkdir -p "$(dirname "$STATUS")" "$ROOT/03_population_structure/logs" "$ROOT/04_snp_indel_gwas/logs"

printf 'status\trunning\nhost\t%s\npid\t%s\nstarted\t%s\n' \
  "$(hostname)" "$$" "$(date '+%F %T %Z')" > "$STATUS"

on_error() {
  local rc=$?
  printf 'status\tfailed\nhost\t%s\npid\t%s\nfailed\t%s\nexit_code\t%s\n' \
    "$(hostname)" "$$" "$(date '+%F %T %Z')" "$rc" > "$STATUS"
  exit "$rc"
}
trap on_error ERR

bash "$ROOT/03_population_structure/scripts/run_population_structure.sh" \
  > "$ROOT/03_population_structure/logs/run_population_structure.log" 2>&1

bash "$ROOT/04_snp_indel_gwas/scripts/run_snp_indel_gwas.sh" \
  > "$ROOT/04_snp_indel_gwas/logs/run_snp_indel_gwas.log" 2>&1

printf 'status\tcompleted\nhost\t%s\npid\t%s\ncompleted\t%s\n' \
  "$(hostname)" "$$" "$(date '+%F %T %Z')" > "$STATUS"
