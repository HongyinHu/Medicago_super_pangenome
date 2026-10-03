#!/usr/bin/env bash
set -euo pipefail
cd path/to/project/37.karyotype_reconstruction/output_Kary/genome_M46_2
trap 'rc=$?; echo RUN_EXIT_STATUS=$rc; exit "$rc"' EXIT
WGDI=path/to/home/anaconda3/envs/wgdi/bin/wgdi
echo STEP=dotplot
"$WGDI" -d total.conf
echo STEP=improved_collinearity
"$WGDI" -icl total.conf
echo STEP=ks
"$WGDI" -ks total.conf
echo STEP=block_information
"$WGDI" -bi total.conf
echo STEP=correspondence
"$WGDI" -c total.conf
echo STEP=block_ks
"$WGDI" -bk total.conf
echo STEP=karyotype_mapping
"$WGDI" -km total.conf
echo STEP=karyotype_plot
"$WGDI" -k total.conf
touch wgdi_complete.flag
echo COMPLETE