#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/N_1.EDTA_single
PAN=$BASE/03_panEDTA
OUTDIR=$BASE/07_LAI_panEDTA_IDfixed_fullpass
LOGROOT=$OUTDIR/logs
PY=path/to/home/anaconda3/envs/EDTA_env/bin/python

mkdir -p "$OUTDIR" "$LOGROOT"

find_intact_passlist() {
  local genome="$1"
  local sample="$2"
  local cand
  for cand in \
    "$PAN/00_inputs/$sample/$genome.mod.EDTA.raw/LTR/$genome.mod.pass.list" \
    "$BASE/02_EDTA_single/$sample/$genome.mod.EDTA.raw/LTR/$genome.mod.pass.list" \
    "$PAN/$genome.mod.EDTA.raw/LTR/$genome.mod.pass.list"; do
    if [ -s "$cand" ]; then
      printf '%s\n' "$cand"
      return 0
    fi
  done
  cand=$(find -L "$BASE/02_EDTA_single/$sample" "$PAN/00_inputs/$sample" \
    \( -path '*/EDTA.raw/LTR/*.pass.list' -o -path '*mod.EDTA.raw/LTR/*.pass.list' \) \
    -type f -size +0c ! -name '*.nmtf.pass.list' 2>/dev/null | sort | head -1)
  if [ -n "$cand" ] && [ -s "$cand" ]; then
    printf '%s\n' "$cand"
    return 0
  fi
  find -L "$BASE/02_EDTA_single/$sample" "$PAN/00_inputs/$sample" \
    \( -path '*/EDTA.raw/LTR/*.pass.list' -o -path '*mod.EDTA.raw/LTR/*.pass.list' \) \
    -type f -size +0c 2>/dev/null | sort | head -1
}

if [ ! -s "$PAN/genome.list" ]; then
  echo "missing $PAN/genome.list" >&2
  exit 2
fi

awk '{print $1}' "$PAN/genome.list" | sed 's#.*/##' > "$OUTDIR/genome.queue"
printf 'sample\toriginal_genome\tmod_genome\tpass_list\tpanEDTA_out\tidfixed_out\tseq_count\tmapped_records\tunchanged_records\tunmapped_records\tstatus\n' > "$OUTDIR/prepare_manifest.tsv.tmp"

while IFS= read -r genome; do
  [ -n "$genome" ] || continue
  sample=${genome%.fa}
  work="$OUTDIR/$sample"
  mkdir -p "$work"

  orig="$PAN/00_inputs/$sample/$genome"
  if [ ! -s "$orig" ]; then
    orig="$PAN/$genome"
  fi
  mod="$PAN/$genome.mod"
  panout="$PAN/$genome.mod.panEDTA.out"
  pass=$(find_intact_passlist "$genome" "$sample" || true)

  if [ ! -s "$orig" ] || [ ! -s "$mod" ] || [ ! -s "$panout" ] || [ -z "${pass:-}" ] || [ ! -s "$pass" ]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\tmissing_input\n' \
      "$sample" "$orig" "$mod" "${pass:-NA}" "$panout" "$work/$genome.mod.panEDTA.IDfixed.out" >> "$OUTDIR/prepare_manifest.tsv.tmp"
    continue
  fi

  ln -sf "$orig" "$work/$genome.mod"
  ln -sf "$pass" "$work/$genome.mod.pass.list"

  "$PY" - "$orig" "$mod" "$panout" "$work/$genome.mod.panEDTA.IDfixed.out" "$work/id_map.mod_to_original.tsv" "$work/rewrite.stats.tsv" <<'PY'
import sys
from pathlib import Path

orig_fa, mod_fa, panout, out_path, map_path, stats_path = map(Path, sys.argv[1:])

def headers(path):
    hs = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                hs.append(line[1:].strip().split()[0])
    return hs

orig = headers(orig_fa)
mod = headers(mod_fa)
if not orig or not mod:
    raise SystemExit(f"empty headers: original={len(orig)} mod={len(mod)}")
if len(orig) != len(mod):
    raise SystemExit(f"header count mismatch: original={len(orig)} mod={len(mod)}")

mapping = dict(zip(mod, orig))
with map_path.open("w", encoding="utf-8") as out:
    out.write("mod_id\toriginal_id\n")
    for m, o in zip(mod, orig):
        out.write(f"{m}\t{o}\n")

total = mapped = unchanged = unmapped = 0
unmapped_ids = {}
with panout.open("r", encoding="utf-8", errors="replace") as inp, out_path.open("w", encoding="utf-8") as out:
    for line in inp:
        stripped = line.strip()
        if not stripped:
            out.write(line)
            continue
        parts = stripped.split()
        if len(parts) >= 5 and parts[0].isdigit():
            total += 1
            q = parts[4]
            if q in mapping:
                parts[4] = mapping[q]
                mapped += 1
            else:
                unchanged += 1
                if q.startswith("_J"):
                    unmapped += 1
                    unmapped_ids[q] = unmapped_ids.get(q, 0) + 1
            out.write("\t".join(parts) + "\n")
        else:
            out.write(line)

with stats_path.open("w", encoding="utf-8") as out:
    out.write("seq_count\tmapped_records\tunchanged_records\tunmapped_records\tunmapped_ids\n")
    out.write(f"{len(mapping)}\t{mapped}\t{unchanged}\t{unmapped}\t{','.join(list(unmapped_ids)[:20])}\n")
PY

  read -r seq_count mapped unchanged unmapped _ < <(awk 'NR==2{print $1, $2, $3, $4, $5}' "$work/rewrite.stats.tsv")
  status=ok
  if [ "${unmapped:-0}" -ne 0 ]; then
    status=unmapped_ids
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$(readlink -f "$orig")" "$(readlink -f "$mod")" "$(readlink -f "$pass")" "$(readlink -f "$panout")" "$work/$genome.mod.panEDTA.IDfixed.out" \
    "${seq_count:-0}" "${mapped:-0}" "${unchanged:-0}" "${unmapped:-0}" "$status" >> "$OUTDIR/prepare_manifest.tsv.tmp"
done < "$OUTDIR/genome.queue"

mv "$OUTDIR/prepare_manifest.tsv.tmp" "$OUTDIR/prepare_manifest.tsv"
echo "Prepared full-pass ID-fixed LAI inputs: $OUTDIR/prepare_manifest.tsv"
