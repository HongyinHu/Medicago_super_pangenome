#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
MANIFEST=download_jing_ye_illumina_RNA_manifest.tsv
LOG=download_jing_ye_illumina_RNA.log
DONE=download_jing_ye_illumina_RNA.done
LOCK=download_jing_ye_illumina_RNA.lock
N_DOWNLOADS=${N_DOWNLOADS:-3}
exec 9>"$LOCK"
if ! flock -n 9; then
  echo "[$(date '+%F %T')] Another download process is running" | tee -a "$LOG"
  exit 1
fi
rm -f "$DONE"
: > "$LOG"
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
download_one(){
  local sample="$1" run="$2" read="$3" url="$4" md5="$5" bytes="$6" final="$7"
  local tmp="${final}.part" got actual
  {
    echo "[$(date '+%F %T')] START ${sample} ${run} read${read} -> ${final}"
    if [ -s "$final" ]; then
      got=$(md5sum "$final" | awk '{print $1}')
      if [ "$got" = "$md5" ]; then
        echo "[$(date '+%F %T')] OK existing ${final}"
        exit 0
      fi
      echo "[$(date '+%F %T')] Existing ${final} md5 mismatch: ${got} != ${md5}; moving aside"
      mv -f "$final" "${final}.bad.$(date +%Y%m%d%H%M%S)"
    fi
    wget -c -nv --tries=0 --timeout=60 --waitretry=20 -O "$tmp" "$url"
    actual=$(stat -c%s "$tmp")
    echo "[$(date '+%F %T')] ${tmp} bytes=${actual} expected=${bytes}"
    if [ "$actual" != "$bytes" ]; then
      echo "ERROR byte_size_mismatch ${tmp}" >&2
      exit 2
    fi
    got=$(md5sum "$tmp" | awk '{print $1}')
    echo "[$(date '+%F %T')] ${tmp} md5=${got} expected=${md5}"
    if [ "$got" != "$md5" ]; then
      echo "ERROR md5_mismatch ${tmp}" >&2
      exit 3
    fi
    mv -f "$tmp" "$final"
    echo "[$(date '+%F %T')] DONE ${final}"
  } >> "$LOG" 2>&1
}
log "START all downloads; N_DOWNLOADS=${N_DOWNLOADS}; host=$(hostname); pwd=$(pwd)"
fail=0
exec 3< <(tail -n +2 "$MANIFEST")
while IFS=$'\t' read -r sample run read url md5 bytes final <&3; do
  download_one "$sample" "$run" "$read" "$url" "$md5" "$bytes" "$final" &
  while [ "$(jobs -rp | wc -l)" -ge "$N_DOWNLOADS" ]; do
    if ! wait -n; then fail=1; fi
  done
done
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
  if ! wait -n; then fail=1; fi
done
if [ "$fail" != 0 ]; then
  log "ERROR at least one download failed"
  exit 10
fi
log "Running gzip integrity test"
gzip -t jingA_1.fastq.gz jingA_2.fastq.gz jingB_1.fastq.gz jingB_2.fastq.gz jingC_1.fastq.gz jingC_2.fastq.gz yeA_1.fastq.gz yeA_2.fastq.gz yeB_1.fastq.gz yeB_2.fastq.gz yeC_1.fastq.gz yeC_2.fastq.gz
md5sum jingA_1.fastq.gz jingA_2.fastq.gz jingB_1.fastq.gz jingB_2.fastq.gz jingC_1.fastq.gz jingC_2.fastq.gz yeA_1.fastq.gz yeA_2.fastq.gz yeB_1.fastq.gz yeB_2.fastq.gz yeC_1.fastq.gz yeC_2.fastq.gz > jing_ye_illumina_RNA.fastq.gz.md5
ls -lh jingA_1.fastq.gz jingA_2.fastq.gz jingB_1.fastq.gz jingB_2.fastq.gz jingC_1.fastq.gz jingC_2.fastq.gz yeA_1.fastq.gz yeA_2.fastq.gz yeB_1.fastq.gz yeB_2.fastq.gz yeC_1.fastq.gz yeC_2.fastq.gz download_jing_ye_illumina_RNA_manifest.tsv jing_ye_illumina_RNA.fastq.gz.md5 >> "$LOG" 2>&1
date '+%F %T' > "$DONE"
log "DONE all downloads"
