#!/usr/bin/env bash
set -euo pipefail
BASE=path/to/project/N_1.EDTA_single
PAN="$BASE/03_panEDTA"
OUT="$BASE/08_TEsorter_EDTA_final"
mkdir -p "$OUT/results"
manifest="$OUT/manifest.tsv"
tmp="$manifest.tmp"
printf 'sample\ttype\tinput_path\tinput_size\tinput_mtime\tout_dir\n' > "$tmp"
add_one(){
  sample="$1"; type="$2"; input="$3"; outdir="$OUT/results/$sample"
  [ -s "$input" ] || { echo "WARN missing input: $input" >&2; return 0; }
  mkdir -p "$outdir"
  size=$(stat -c %s "$input")
  mtime=$(stat -c '%y' "$input" | cut -d. -f1)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$type" "$input" "$size" "$mtime" "$outdir" >> "$tmp"
}
add_one panEDTA_combined panEDTA_final "$PAN/genome.list.panEDTA.TElib.fa"
find "$PAN" -maxdepth 1 -type f -name '*.fa.mod.EDTA.TElib.fa' | sort | while read -r f; do
  b=$(basename "$f")
  sample=${b%.fa.mod.EDTA.TElib.fa}
  add_one "$sample" per_genome_panEDTA_final "$f"
done
mv "$tmp" "$manifest"
echo "manifest=$manifest"
echo "tasks=$(( $(wc -l < "$manifest") - 1 ))"
