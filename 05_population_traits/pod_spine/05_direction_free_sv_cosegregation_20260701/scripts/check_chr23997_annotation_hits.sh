#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project
OUT="$ROOT/16.T2T_ref_function_anno/output"
cd "$OUT"

echo "CANDIDATE_FILES"
find . -maxdepth 5 -type f | grep -Ei '(tair|blast|swiss|interpro|anno|tsv|txt|xls|m8|out)' | sed -n '1,240p'

echo "CHR23997_ALL_HITS"
find . -maxdepth 6 -type f | xargs -r grep -H -i -E 'Chr23997|Chr239971|Chr23997\.1' 2>/dev/null | sed -n '1,240p'

echo "INTERPRO_MAPPING"
find . -maxdepth 6 -type f | xargs -r grep -H -i -E 'IPR004333|IPR036893|Squamosa|SBP|SPL' 2>/dev/null | sed -n '1,240p'
