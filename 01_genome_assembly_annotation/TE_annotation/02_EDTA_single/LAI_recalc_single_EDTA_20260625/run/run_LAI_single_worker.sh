#!/usr/bin/env bash
set -uo pipefail
NODE=${1:-$(hostname -s)}
THREADS=${2:-30}
OUT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
MANIFEST="$OUT/manifest_single_EDTA_LAI.tsv"
LOGDIR="$OUT/logs/$NODE"
mkdir -p "$LOGDIR" "$OUT/work" "$OUT/status"
DRIVER="$LOGDIR/driver.log"
{
  echo "[$(date '+%F %T')] worker_start node=$NODE host=$(hostname) threads=$THREADS out=$OUT"
  source path/to/home/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || true
  conda activate EDTA_env
  echo "[$(date '+%F %T')] LAI=$(command -v LAI || true)"
} >> "$DRIVER" 2>&1

run_one() {
  local sample="$1" fa="$2" pass="$3" rmout="$4"
  local lock="$OUT/status/$sample.lock"
  local done="$OUT/status/$sample.done"
  local failed="$OUT/status/$sample.failed"
  [ -s "$done" ] && return 0
  [ -s "$failed" ] && return 0
  if ! mkdir "$lock" 2>/dev/null; then
    return 0
  fi
  printf '%s\t%s\t%s\n' "$NODE" "$(hostname)" "$(date '+%F %T')" > "$lock/info.tsv"
  local wd="$OUT/work/$sample"
  local log="$LOGDIR/$sample.log"
  mkdir -p "$wd"
  {
    echo "[$(date '+%F %T')] sample_start=$sample node=$NODE host=$(hostname) threads=$THREADS"
    echo "fa=$fa"
    echo "pass_list=$pass"
    echo "rmout=$rmout"
    cd "$wd" || exit 2
    ln -sfn "$fa" genome.fa.mod
    ln -sfn "$pass" intact.pass.list
    ln -sfn "$rmout" all.mod.out
    rm -f all.mod.out.LAI all.mod.out.LAI.* genome.fa.mod.LAI* intact.pass.list.LAI* 2>/dev/null || true
    if command -v /usr/bin/time >/dev/null 2>&1; then
      /usr/bin/time -v LAI -genome genome.fa.mod -intact intact.pass.list -all all.mod.out -t "$THREADS"
    else
      LAI -genome genome.fa.mod -intact intact.pass.list -all all.mod.out -t "$THREADS"
    fi
    rc=$?
    echo "[$(date '+%F %T')] sample_end=$sample rc=$rc"
    if [ "$rc" -eq 0 ] && [ -s all.mod.out.LAI ]; then
      printf 'sample\tnode\thost\tthreads\tfinished\tworkdir\tlai_file\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$NODE" "$(hostname)" "$THREADS" "$(date '+%F %T')" "$wd" "$wd/all.mod.out.LAI" > "$done"
      rm -f "$failed"
    else
      printf 'sample\tnode\thost\tthreads\tfailed_time\trc\tworkdir\tlog\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$NODE" "$(hostname)" "$THREADS" "$(date '+%F %T')" "$rc" "$wd" "$log" > "$failed"
    fi
  } > "$log" 2>&1
  rm -rf "$lock"
}

while IFS=$'\t' read -r sample fa pass rmout status; do
  [ "$sample" = "sample" ] && continue
  [ "$status" = "ready" ] || continue
  run_one "$sample" "$fa" "$pass" "$rmout"
done < "$MANIFEST"

echo "[$(date '+%F %T')] worker_done node=$NODE host=$(hostname)" >> "$DRIVER"
