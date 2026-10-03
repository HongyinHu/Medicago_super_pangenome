#!/usr/bin/env bash
set -euo pipefail
ROOT=path/to/project
RUN="$ROOT/N_5.pod_coil/01_SHP_homology_and_coil_cosegregation"
ENV=path/to/home/anaconda3/envs/biosofeware/bin
cd "$RUN/trees"
"$ENV/iqtree2" -s "$RUN/alignments/SHP_reference_tree.pep.aln.fa" -m MFP -B 1000 -alrt 1000 \
  -T 4 --prefix SHP_reference_fixed4 >"$RUN/logs/SHP_reference_fixed4.iqtree.log" 2>&1
touch "$RUN/SHP_reference_tree.done"
