#!/usr/bin/env bash
set -euo pipefail
NODE_LABEL=${1:-$(hostname)}
THREADS=${THREADS:-12}
BASE=path/to/project/N_1.EDTA_single
OUT="$BASE/08_TEsorter_EDTA_final"
MANIFEST="$OUT/manifest.tsv"
LOGDIR="$OUT/logs/$NODE_LABEL"
mkdir -p "$LOGDIR" "$OUT/task_locks"
log="$LOGDIR/driver.log"
{
  echo "[$(date)] START worker host=$(hostname) node_label=$NODE_LABEL threads=$THREADS"
  source ~/.bashrc 2>/dev/null || true
  conda activate EDTA_env
  echo "[$(date)] TEsorter=$(command -v TEsorter)"
  tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample type input size mtime outdir; do
    [ -n "$sample" ] || continue
    done_marker="$outdir/.done"
    fail_marker="$outdir/.failed"
    lockdir="$OUT/task_locks/$sample.lock"
    if [ -s "$done_marker" ]; then
      echo "[$(date)] SKIP done $sample"
      continue
    fi
    if ! mkdir "$lockdir" 2>/dev/null; then
      echo "[$(date)] SKIP locked $sample"
      continue
    fi
    trap 'rm -rf "$lockdir"' EXIT
    {
      echo "host=$(hostname)"
      echo "node_label=$NODE_LABEL"
      echo "pid=$$"
      echo "start=$(date '+%F %T %Z')"
      echo "input=$input"
    } > "$lockdir/info.txt"
    mkdir -p "$outdir"
    rm -f "$fail_marker"
    cd "$outdir"
    ln -sf "$input" input.TElib.fa
    echo "[$(date)] RUN $sample type=$type input_size=$size"
    set +e
    /usr/bin/time -v TEsorter input.TElib.fa -db rexdb-plant -p "$THREADS" -pre "$sample.rexdb-plant" -tmp "$outdir/tmp" > "$outdir/TEsorter.stdout.log" 2> "$outdir/TEsorter.stderr.log"
    rc=$?
    set -e
    if [ "$rc" -eq 0 ]; then
      cls_count=$(find "$outdir" -maxdepth 1 -type f \( -name '*.cls.tsv' -o -name '*.cls.pep' -o -name '*.domtbl' -o -name '*.gff3' \) | wc -l)
      {
        echo "status=ok"
        echo "sample=$sample"
        echo "type=$type"
        echo "input=$input"
        echo "threads=$THREADS"
        echo "host=$(hostname)"
        echo "finished=$(date '+%F %T %Z')"
        echo "output_files=$cls_count"
        find "$outdir" -maxdepth 1 -type f -printf '%f\t%s\n' | sort
      } > "$done_marker"
      echo "[$(date)] DONE $sample outputs=$cls_count"
    else
      {
        echo "status=failed"
        echo "sample=$sample"
        echo "type=$type"
        echo "input=$input"
        echo "threads=$THREADS"
        echo "host=$(hostname)"
        echo "failed=$(date '+%F %T %Z')"
        echo "exit_code=$rc"
        tail -80 "$outdir/TEsorter.stderr.log" 2>/dev/null || true
      } > "$fail_marker"
      echo "[$(date)] FAILED $sample rc=$rc"
    fi
    rm -rf "$lockdir"
    trap - EXIT
    cd "$OUT"
  done
  echo "[$(date)] END worker host=$(hostname) node_label=$NODE_LABEL"
} >> "$log" 2>&1
