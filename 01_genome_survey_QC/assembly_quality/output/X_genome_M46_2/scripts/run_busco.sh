#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/6.genome_quality_assess/output/X_genome_M46_2
export PATH=path/to/home/anaconda3/envs/busco_env/bin:$PATH
BUSCO_PY=path/to/home/anaconda3/envs/busco_env/bin/python
BUSCO=path/to/home/anaconda3/envs/busco_env/bin/busco
LINEAGE=path/to/home/busco_version/embryophyta_1614_odb10
cd "$BASE/assembly_evaluate/BUSCO"
"$BUSCO_PY" "$BUSCO" -i "$BASE/data/M46_genome_Chr_reorder.genome.ctg.fa" -l "$LINEAGE" -o busco_assess --out_path . -c 16 -m genome --offline -f > busco.log 2>&1
cd "$BASE/annotation_evaluate/BUSCO"
"$BUSCO_PY" "$BUSCO" -i "$BASE/data/finally.genome.anno.pep" -l "$LINEAGE" -o busco_pep --out_path . -c 16 -m proteins --offline -f > busco.log 2>&1
