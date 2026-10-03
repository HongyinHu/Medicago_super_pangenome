#!/usr/bin/env bash
set -euo pipefail
OUT=path/to/project/N_2.DNA_survy/output
PID_FILE="$OUT/03_estimates/kmer_survey.pid"
echo "== PID =="
if [ -f "$PID_FILE" ]; then
  PID=$(cat "$PID_FILE")
  echo "$PID"
  ps -p "$PID" -o pid,ppid,stat,etime,cmd || true
else
  echo "No PID file"
fi

echo "\n== Related processes =="
pgrep -af '01_run_jellyfish_kmer_survey|jellyfish|gzip -cd' || true

echo "\n== Sample counts =="
TOTAL=$(tail -n +2 "$OUT/01_sample_table/samples.tsv" | wc -l)
DONE=$(find "$OUT/02_jellyfish/k21" -maxdepth 1 -name '*.histo' 2>/dev/null | wc -l)
echo "done_k21=$DONE / total=$TOTAL"

echo "\n== Last status lines =="
[ -f "$OUT/03_estimates/pipeline_status.tsv" ] && tail -20 "$OUT/03_estimates/pipeline_status.tsv" || true

echo "\n== Estimate table =="
if [ -f "$OUT/03_estimates/genome_size_estimates.tsv" ]; then
  column -t -s $'\t' "$OUT/03_estimates/genome_size_estimates.tsv" | sed -n '1,25p'
else
  echo "No estimates yet"
fi

echo "\n== Latest sample log =="
latest=$(ls -t "$OUT"/logs/*.k21.log 2>/dev/null | head -1 || true)
[ -n "${latest:-}" ] && { echo "$latest"; tail -30 "$latest"; } || true
