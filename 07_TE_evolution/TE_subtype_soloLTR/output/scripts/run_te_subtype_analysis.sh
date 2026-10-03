#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-path/to/project/39.TE_type_soloLTR}"
DATA_DIR="${DATA_DIR:-$BASE/data}"
OUT_DIR="${OUT_DIR:-$BASE/output}"
THREADS="${THREADS:-16}"
DB="${DB:-rexdb-plant}"
CONDA_SH="${CONDA_SH:-path/to/home/anaconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-EDTA_env}"

mkdir -p "$OUT_DIR/tesorter" "$OUT_DIR/logs"

source "$CONDA_SH"
conda activate "$CONDA_ENV"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
summary_script="$script_dir/summarize_te_subtypes.py"

if [[ ! -s "$summary_script" ]]; then
  echo "Missing summary script: $summary_script" >&2
  exit 1
fi

echo "[$(date '+%F %T')] host=$(hostname) env=$CONDA_ENV db=$DB threads=$THREADS"

for edta_dir in "$DATA_DIR"/*_EDTA; do
  [[ -e "$edta_dir" ]] || continue
  species="$(basename "$edta_dir")"
  species="${species%_EDTA}"
  telib="$(find -L "$edta_dir" -maxdepth 1 -type f -name '*.EDTA.TElib.fa' | head -n 1)"

  if [[ -z "$telib" ]]; then
    echo "[$(date '+%F %T')] WARN no TElib found under $edta_dir" >&2
    continue
  fi

  species_out="$OUT_DIR/tesorter/$species"
  mkdir -p "$species_out"
  prefix="$species_out/$species.$DB"
  cls="$prefix.cls.tsv"

  if [[ -s "$cls" ]]; then
    echo "[$(date '+%F %T')] skip existing $cls"
  else
    echo "[$(date '+%F %T')] TEsorter $species"
    rm -rf "$species_out/tmp"
    TEsorter "$telib" \
      -db "$DB" \
      -p "$THREADS" \
      -tmp "$species_out/tmp" \
      -pre "$prefix" \
      > "$species_out/TEsorter.log" 2>&1
  fi
done

echo "[$(date '+%F %T')] summarize"
python "$summary_script" \
  --base "$BASE" \
  --db "$DB" \
  --output-dir "$OUT_DIR"

echo "[$(date '+%F %T')] done"
