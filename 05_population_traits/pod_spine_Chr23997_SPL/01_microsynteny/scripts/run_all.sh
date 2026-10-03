#!/usr/bin/env bash
set -euo pipefail

ROOT="path/to/project/41.Chr23997_realName"
DATA="$ROOT/00_data"
OUT="$ROOT/01_microsynteny"
SCRIPT_DIR="$OUT/scripts"
LOG="$OUT/commands_log.txt"

mkdir -p "$OUT"/{blastdb,blast,logs,scripts}

{
  echo "# Chr23997 microsynteny analysis"
  echo "date: $(date -Iseconds)"
  echo "host: $(hostname)"
  echo "root: $ROOT"
  echo "data: $DATA"
  echo
  echo "## Tool discovery"
} > "$LOG"

find_tool() {
  local exe="$1"
  if command -v "$exe" >/dev/null 2>&1; then
    command -v "$exe"
    return 0
  fi
  for d in \
    path/to/home/anaconda3/bin \
    path/to/home/anaconda3/envs/*/bin \
    path/to/home/bin \
    path/to/home/sofeware/*/bin \
    path/to/home/sofeware/*; do
    if [ -x "$d/$exe" ]; then
      echo "$d/$exe"
      return 0
    fi
  done
  return 1
}

BLASTP="$(find_tool blastp || true)"
MAKEBLASTDB="$(find_tool makeblastdb || true)"
PYTHON_BIN="$(find_tool python3 || command -v python3)"
RSCRIPT_BIN="$(find_tool Rscript || command -v Rscript)"

{
  echo "blastp: ${BLASTP:-NOT_FOUND}"
  echo "makeblastdb: ${MAKEBLASTDB:-NOT_FOUND}"
  echo "python3: $PYTHON_BIN"
  echo "Rscript: $RSCRIPT_BIN"
  echo
  echo "## Commands"
  printf 'BLASTP=%q MAKEBLASTDB=%q %q %q --data %q --out %q\n' \
    "$BLASTP" "$MAKEBLASTDB" "$PYTHON_BIN" "$SCRIPT_DIR/run_microsynteny.py" "$DATA" "$OUT"
} >> "$LOG"

if [ -z "$BLASTP" ] || [ -z "$MAKEBLASTDB" ]; then
  echo "ERROR: blastp or makeblastdb was not found on this node." | tee -a "$LOG"
  exit 2
fi

BLASTP="$BLASTP" MAKEBLASTDB="$MAKEBLASTDB" \
  "$PYTHON_BIN" "$SCRIPT_DIR/run_microsynteny.py" --data "$DATA" --out "$OUT" \
  2>&1 | tee "$OUT/logs/run_microsynteny.log"

{
  printf '%q %q %q %q\n' "$RSCRIPT_BIN" "$SCRIPT_DIR/plot_microsynteny.R" "$OUT" "$OUT/microsynteny_plot"
} >> "$LOG"

"$RSCRIPT_BIN" "$SCRIPT_DIR/plot_microsynteny.R" "$OUT" "$OUT/microsynteny_plot" \
  2>&1 | tee "$OUT/logs/plot_microsynteny.log"

{
  echo
  echo "## Outputs"
  ls -lh "$OUT"/neighboring_gene_table.tsv \
    "$OUT"/reciprocal_blast_results.tsv \
    "$OUT"/microsynteny_plot.pdf \
    "$OUT"/microsynteny_plot.png \
    "$OUT"/interpretation_report.md
} >> "$LOG"

