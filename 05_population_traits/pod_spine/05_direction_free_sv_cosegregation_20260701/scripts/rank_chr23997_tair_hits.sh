#!/usr/bin/env bash
set -euo pipefail

TAB=path/to/project/16.T2T_ref_function_anno/output/blast_tair/T2T_tait/ref_Msa.T2T_ctg.pep.tab
ANNO=path/to/project/16.T2T_ref_function_anno/output/blast_tair/T2T_tait/ref_Msa.T2T_ctg.pep.tab.anno

echo "BEST_TAIR_HITS_BY_BITSCORE"
awk '$2=="Chr239971"{print $0}' "$TAB" | sort -k12,12gr | head -20

echo "ANNOTATED_LINES_FOR_TOP_SPL"
grep -E 'Chr239971.*(SPL2|SPL10|SPL11|SPL8|SPL14)' "$ANNO" | sed -n '1,80p'

echo "INTERPRO_DOMAIN"
grep -w Chr23997 path/to/project/16.T2T_ref_function_anno/output/interproscan/ref_Msa.T2T_ctg.pep.domains.txt || true
