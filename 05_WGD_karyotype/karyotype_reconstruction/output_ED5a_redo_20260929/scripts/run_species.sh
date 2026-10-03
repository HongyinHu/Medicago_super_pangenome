#!/usr/bin/env bash
# Run WGDI -c -> one-to-one -> -km -> -k for one staged species directory.
set -euo pipefail
OUT=path/to/project/37.karyotype_reconstruction/output_ED5a_redo_20260929
ENV=path/to/home/anaconda3/envs/wgdi/bin
RUN=${RUN:-02_wgdi}
lab=$1
cd "$OUT/$RUN/$lab"
rm -f status.done status.failed
trap 'rc=$?; [ $rc -ne 0 ] && touch status.failed; exit $rc' EXIT
{
  echo "START $(date '+%F %T') host=$(hostname)"
  "$ENV/wgdi" -c c.conf
  "$ENV/python" "$OUT/scripts/one_to_one.py" "${lab}_aak.ortholog_blocks.csv" "${lab}_aak.ortholog_blocks.1to1.csv"
  "$ENV/wgdi" -km km.conf
  "$ENV/wgdi" -k k.conf
  echo "END $(date '+%F %T')"
} > run.log 2>&1
touch status.done
