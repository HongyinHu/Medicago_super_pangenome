#!/usr/bin/env bash
set -euo pipefail
sample="$1"
threads="${THREADS:-48}"
BASE=path/to/project/N_1.EDTA_single
FDIR=$BASE/07_LAI_panEDTA_IDfixed_fullpass
[ -d $BASE/04_LAI_panEDTA_IDfixed_fullpass ] && [ ! -L $BASE/04_LAI_panEDTA_IDfixed_fullpass ] && FDIR=$BASE/04_LAI_panEDTA_IDfixed_fullpass
FAST=$FDIR/fast_rerun_20260630_474ab
LOG=$FDIR/logs/fast_rerun_20260630_474ab
orig="$FDIR/$sample"
work="$FAST/$sample"
mkdir -p "$work" "$LOG"
cd "$work"
for f in "$sample.fa.mod" "$sample.fa.mod.pass.list" "$sample.fa.mod.panEDTA.IDfixed.out"; do
  [ -e "$orig/$f" ] || { echo "missing $orig/$f" >&2; exit 2; }
  [ -e "$f" ] || ln -s "$orig/$f" "$f"
done
export PATH=path/to/home/anaconda3/envs/EDTA_env/bin:path/to/home/anaconda3/envs/EDTA_env/share/EDTA:path/to/home/anaconda3/envs/EDTA_env/share/LTR_retriever:$PATH
export PYTHONNOUSERSITE=1
start=$(date '+%F %T')
echo "[$start] FAST_LAI_START sample=$sample host=$(hostname) threads=$threads work=$work" | tee -a "$LOG/driver.log"
LAI -genome "$sample.fa.mod" -intact "$sample.fa.mod.pass.list" -all "$sample.fa.mod.panEDTA.IDfixed.out" -t "$threads" > "$LOG/$sample.stdout.$(date +%Y%m%d_%H%M%S).log" 2> "$LOG/$sample.stderr.$(date +%Y%m%d_%H%M%S).log"
if [ -s "$work/$sample.fa.mod.panEDTA.IDfixed.out.LAI" ]; then
  cp -f "$work/$sample.fa.mod.panEDTA.IDfixed.out.LAI" "$orig/$sample.fa.mod.panEDTA.IDfixed.out.LAI.fastcopy.tmp"
  mv -f "$orig/$sample.fa.mod.panEDTA.IDfixed.out.LAI.fastcopy.tmp" "$orig/$sample.fa.mod.panEDTA.IDfixed.out.LAI"
  echo "[$(date '+%F %T')] FAST_LAI_DONE sample=$sample copied_to=$orig/$sample.fa.mod.panEDTA.IDfixed.out.LAI" | tee -a "$LOG/driver.log"
else
  echo "[$(date '+%F %T')] FAST_LAI_NO_FINAL sample=$sample" | tee -a "$LOG/driver.log"
  exit 1
fi
