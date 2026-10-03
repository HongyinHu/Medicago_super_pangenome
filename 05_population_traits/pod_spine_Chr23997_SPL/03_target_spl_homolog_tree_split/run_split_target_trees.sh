#!/usr/bin/env bash
set -e
ROOT=path/to/project/41.Chr23997_realName
OUT=$ROOT/03_target_spl_homolog_tree_split
python3 $ROOT/03_target_spl_homolog_tree_split/scripts/build_split_target_spl_trees.py \
  --data $ROOT/00_data \
  --summary $ROOT/02_target_spl_homolog_tree/target_spl_homolog_summary.tsv \
  --target $ROOT/00_data/Target.SPL.pep.fa \
  --out $OUT \
  --threads 8
