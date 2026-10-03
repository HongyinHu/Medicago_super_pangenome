#!/usr/bin/env bash
set -euo pipefail

ROOT=path/to/project
cd "$ROOT"

echo "PWD"
pwd

echo "N4_00_DATA_FILES"
find N_4.pod_spiny/00_data -maxdepth 4 -type f | sed -n '1,160p'

echo "FUNCTION_ANNOTATION_DIRS"
find . -maxdepth 3 -type d \( -name '*function*' -o -name '*anno*' -o -name '*annotation*' \) | sort | sed -n '1,120p'

echo "CHR23997_GFF_HITS"
find N_4.pod_spiny/00_data -type f \( -name '*.gff' -o -name '*.gff3' -o -name '*.gtf' \) -print0 \
  | xargs -0 -r grep -H -w 'Chr23997' \
  | sed -n '1,80p'

echo "CHR23997_ANNOTATION_HITS"
find 16.T2T_ref_function_anno N_4.pod_spiny/00_data -type f 2>/dev/null \
  | grep -Ei '\.(tsv|txt|csv|gff|gff3|faa|pep|fa|anno|out)$' \
  | xargs -r grep -H -i -w 'Chr23997' 2>/dev/null \
  | sed -n '1,160p'

echo "SPL_HITS_NEAR_CANDIDATE_DATA"
find N_4.pod_spiny/00_data -type f 2>/dev/null \
  | xargs -r grep -H -i -E 'Chr23997|SPL|SBP|SQUAMOSA' 2>/dev/null \
  | sed -n '1,160p'
