#!/usr/bin/env bash
set -euo pipefail
WORK=path/to/project/N_3.call_SV/09_strict_ortholog_expression_SV_distance_20260705
ALL_CDS=path/to/project/1.orthology_family_2/data/pangenome_data/CDS_dir/all.sample.cds
SALMON=path/to/home/anaconda3/envs/biosofeware/bin/salmon
PYTHON=path/to/home/anaconda3/envs/biosofeware/bin/python
RSCRIPT=path/to/home/anaconda3/envs/biosofeware/bin/Rscript
THREADS_INDEX=${THREADS_INDEX:-8}
THREADS_QUANT=${THREADS_QUANT:-6}
JOBS=${JOBS:-3}
cd "$WORK"
mkdir -p reference index quants logs tables figures summary

echo "[$(date '+%F %T')] strict ortholog expression-SV started" | tee logs/run_strict_ortholog_expression_sv.log
awk 'NR>1{print $1"\t"$3}' config/species_map.tsv | while IFS=$'\t' read -r group prefix; do
  out="reference/${group}.cds.fa"
  if [ ! -s "$out" ]; then
    awk -v p="${prefix}|" 'BEGIN{keep=0} /^>/{h=substr($0,2); keep=(index(h,p)==1)} keep{print}' "$ALL_CDS" > "$out"
  fi
  n=$(grep -c '^>' "$out" || true)
  echo -e "$group\t$prefix\t$n" >> logs/reference_cds_counts.tmp
  if [ "$n" -eq 0 ]; then echo "ERROR: no CDS for $group prefix $prefix" >&2; exit 1; fi
done
mv logs/reference_cds_counts.tmp summary/reference_cds_counts.tsv

awk 'NR>1{print $1"\t"$3"\t"$4}' config/samples.tsv | while IFS=$'\t' read -r sample r1 r2; do
  [ -s "$r1" ] || { echo "missing read1 for $sample: $r1" >&2; exit 1; }
  [ -s "$r2" ] || { echo "missing read2 for $sample: $r2" >&2; exit 1; }
done

tail -n +2 config/species_map.tsv | while IFS=$'\t' read -r group ortho prefix pav; do
  idx="index/${group}"
  if [ ! -s "$idx/versionInfo.json" ]; then
    echo "[$(date '+%F %T')] salmon index $group" | tee -a logs/run_strict_ortholog_expression_sv.log
    "$SALMON" index -t "reference/${group}.cds.fa" -i "$idx" -p "$THREADS_INDEX" > "logs/salmon_index.${group}.log" 2>&1
  fi
done

run_one() {
  sample="$1"; group="$2"; r1="$3"; r2="$4"
  out="quants/${sample}"
  if [ -s "$out/quant.sf" ]; then
    echo "[$(date '+%F %T')] skip existing quant $sample"
    return 0
  fi
  rm -rf "$out.tmp"
  echo "[$(date '+%F %T')] salmon quant $sample $group"
  "$SALMON" quant -i "index/${group}" -l A -1 "$r1" -2 "$r2" -p "$THREADS_QUANT" --validateMappings -o "$out.tmp" > "logs/salmon_quant.${sample}.log" 2>&1
  rm -rf "$out"
  mv "$out.tmp" "$out"
}
export -f run_one
export SALMON WORK THREADS_QUANT
awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4}' config/samples.tsv | xargs -P "$JOBS" -n 4 bash -c 'run_one "$@"' _

"$PYTHON" scripts/strict_ortholog_expression_sv.py > logs/strict_ortholog_expression_sv.python.log 2>&1
"$RSCRIPT" scripts/plot_strict_ortholog_expression.R > logs/plot_strict_ortholog_expression.R.log 2>&1

date '+%F %T' > summary/strict_ortholog_expression_sv.done
echo "[$(date '+%F %T')] strict ortholog expression-SV done" | tee -a logs/run_strict_ortholog_expression_sv.log
