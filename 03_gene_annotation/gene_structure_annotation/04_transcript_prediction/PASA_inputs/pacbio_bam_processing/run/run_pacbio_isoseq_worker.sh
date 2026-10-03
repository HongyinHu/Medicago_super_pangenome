#!/usr/bin/env bash
set -uE -o pipefail

BASE=path/to/project/N_1.coding_gene_anno
OUT="$BASE/04_transcript_prediction/PASA_inputs/pacbio_bam_processing"
MAN="$OUT/manifests/pacbio_bam_manifest.tsv"

THREADS=${THREADS:-24}
NODE=$(hostname)

CCS=path/to/home/anaconda3/envs/genome_repeat/bin/ccs
LIMA="$OUT/tools/lima_2.9.0_package/bin/lima"
ISOSEQ=path/to/home/anaconda3/envs/Isoseq3/bin/isoseq
SAMTOOLS=path/to/home/anaconda3/envs/braker3_env/bin/samtools

mkdir -p "$OUT"/{logs,work,locks,status,final}
DRIVER="$OUT/logs/worker_${NODE}_$$_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$DRIVER") 2>&1

echo "[START] $(date '+%F %T %Z') node=$NODE pid=$$ threads=$THREADS"
echo "[TOOLS] $CCS | $LIMA | $ISOSEQ | $SAMTOOLS"

run_duplicate_links() {
  awk -F'\t' 'NR > 1 && $8 == "duplicate_ready" {print $1"\t"$9}' "$MAN" | while IFS=$'\t' read -r sample note; do
    [ -n "$sample" ] || continue
    [ ! -e "$OUT/status/${sample}.done" ] || continue
    source_sample=${note#duplicate_of=}
    [ -n "$source_sample" ] || continue
    [ -e "$OUT/status/${source_sample}.done" ] || continue
    lock="$OUT/locks/${sample}.lock"
    if mkdir "$lock" 2>/dev/null; then
      {
        echo "$NODE $$ $(date '+%F %T %Z')" > "$lock/owner"
        src="$OUT/work/$source_sample"
        dst="$OUT/work/$sample"
        mkdir -p "$dst"
        for suffix in ccs.bam ccs.bam.pbi fl.bam fl.bam.pbi flnc.bam flnc.bam.pbi transcripts.bam transcripts.bam.pbi isoseq_transcripts.fa; do
          if [ -e "$src/${source_sample}.${suffix}" ]; then
            ln -sfn "../$source_sample/${source_sample}.${suffix}" "$dst/${sample}.${suffix}"
          fi
        done
        if [ -e "$src/${source_sample}.isoseq_transcripts.fa" ]; then
          ln -sfn "../work/$sample/${sample}.isoseq_transcripts.fa" "$OUT/final/${sample}.isoseq_transcripts.fa"
          echo "duplicate_of=$source_sample linked_at=$(date '+%F %T %Z')" > "$OUT/status/${sample}.done"
          echo "[DUP_DONE] $sample duplicate_of=$source_sample"
        fi
      } || {
        echo "duplicate link failed at $(date '+%F %T %Z')" > "$OUT/status/${sample}.failed"
      }
      rm -rf "$lock"
    fi
  done
}

process_sample() {
  local sample="$1" bam="$2" primer="$3"
  local lock="$OUT/locks/${sample}.lock"
  [ ! -e "$OUT/status/${sample}.done" ] || return 1
  [ ! -e "$OUT/status/${sample}.failed" ] || return 1
  if ! mkdir "$lock" 2>/dev/null; then
    return 1
  fi

  echo "$NODE $$ $(date '+%F %T %Z')" > "$lock/owner"
  local od="$OUT/work/$sample"
  mkdir -p "$od"
  local sample_log="$OUT/logs/${sample}.${NODE}.$(date +%Y%m%d_%H%M%S).log"
  echo "[CLAIM] $sample bam=$bam primer=$primer log=$sample_log"

  (
    set -euo pipefail
    echo "[SAMPLE_START] $(date '+%F %T %Z') sample=$sample node=$NODE threads=$THREADS"
    echo "[INPUT] bam=$bam"
    echo "[INPUT] primer=$primer"
    local_tmp_root=${LOCAL_TMP_ROOT:-/tmp/${USER}_pacbio_isoseq}
    local_tmp="$local_tmp_root/${sample}.${NODE}.${BASHPID}"
    mkdir -p "$local_tmp"
    cleanup_local_tmp() {
      rm -rf -- "$local_tmp"
    }
    trap cleanup_local_tmp EXIT
    export TMPDIR="$local_tmp"
    export TMP="$local_tmp"
    export TEMP="$local_tmp"
    echo "[LOCAL_TMP] $local_tmp"

    if [ ! -s "$od/${sample}.ccs.bam" ]; then
      "$CCS" -j "$THREADS" --min-rq 0.9 \
        --report-file "$od/${sample}.ccs.report.txt" \
        --report-json "$od/${sample}.ccs.report.json" \
        "$bam" "$od/${sample}.ccs.bam"
    fi

    fl_input="$od/${sample}.fl.bam"
    if [ ! -s "$fl_input" ]; then
      fl_input=$(find "$od" -maxdepth 1 -type f -name "${sample}.fl.*--*.bam" -size +0c | sort | head -1 || true)
    fi
    if [ -z "${fl_input:-}" ] || [ ! -s "$fl_input" ]; then
      "$LIMA" -j "$THREADS" --isoseq --peek-guess "$od/${sample}.ccs.bam" "$primer" "$od/${sample}.fl.bam"
      fl_input="$od/${sample}.fl.bam"
      if [ ! -s "$fl_input" ]; then
        fl_input=$(find "$od" -maxdepth 1 -type f -name "${sample}.fl.*--*.bam" -size +0c | sort | head -1 || true)
      fi
    fi
    if { [ -z "${fl_input:-}" ] || [ ! -s "$fl_input" ]; } && [ -s "$od/${sample}.fl.consensusreadset.xml" ]; then
      fl_input="$od/${sample}.fl.consensusreadset.xml"
    fi
    if [ -z "${fl_input:-}" ] || [ ! -s "$fl_input" ]; then
      echo "[ERROR] missing FL Iso-Seq input after lima for $sample"
      exit 10
    fi
    echo "[FL_INPUT] $fl_input"

    if [ ! -s "$od/${sample}.flnc.bam" ]; then
      "$ISOSEQ" refine -j "$THREADS" --require-polya "$fl_input" "$primer" "$od/${sample}.flnc.bam"
    fi
    if [ ! -s "$od/${sample}.flnc.bam" ]; then
      echo "[ERROR] missing FLNC BAM for $sample"
      exit 11
    fi

    if [ ! -s "$od/${sample}.transcripts.bam" ]; then
      "$ISOSEQ" cluster2 -j "$THREADS" "$od/${sample}.flnc.bam" "$od/${sample}.transcripts.bam"
    fi
    if [ ! -s "$od/${sample}.transcripts.bam" ]; then
      echo "[ERROR] missing clustered transcripts BAM for $sample"
      exit 12
    fi

    if [ ! -s "$od/${sample}.isoseq_transcripts.fa" ]; then
      "$SAMTOOLS" fasta "$od/${sample}.transcripts.bam" > "$od/${sample}.isoseq_transcripts.fa.tmp"
      mv "$od/${sample}.isoseq_transcripts.fa.tmp" "$od/${sample}.isoseq_transcripts.fa"
    fi

    ln -sfn "../work/$sample/${sample}.isoseq_transcripts.fa" "$OUT/final/${sample}.isoseq_transcripts.fa"
    nseq=$(grep -c '^>' "$od/${sample}.isoseq_transcripts.fa" || true)
    bytes=$(stat -Lc '%s' "$od/${sample}.isoseq_transcripts.fa" 2>/dev/null || echo 0)
    if [ "$nseq" -eq 0 ] || [ "$bytes" -eq 0 ]; then
      echo "[ERROR] empty Iso-Seq transcript FASTA for $sample"
      exit 13
    fi
    {
      echo "sample=$sample"
      echo "node=$NODE"
      echo "finished=$(date '+%F %T %Z')"
      echo "threads=$THREADS"
      echo "transcripts=$nseq"
      echo "fasta_bytes=$bytes"
    } > "$OUT/status/${sample}.done"
    echo "[SAMPLE_DONE] $(date '+%F %T %Z') sample=$sample transcripts=$nseq fasta_bytes=$bytes"
  ) > "$sample_log" 2>&1
  rc=$?

  if [ "$rc" -ne 0 ]; then
    {
      echo "sample=$sample"
      echo "node=$NODE"
      echo "failed=$(date '+%F %T %Z')"
      echo "rc=$rc"
      echo "log=$sample_log"
    } > "$OUT/status/${sample}.failed"
    echo "[SAMPLE_FAILED] sample=$sample rc=$rc log=$sample_log"
  fi

  rm -rf "$lock"
  return 0
}

while true; do
  run_duplicate_links
  claimed=0
  while IFS=$'\t' read -r sample input_bam real_bam size_bytes pbi primer real_primer state note; do
    [ "$sample" != "sample" ] || continue
    [ "$state" = "ready_full" ] || continue
    if process_sample "$sample" "$input_bam" "$primer"; then
      claimed=1
      break
    fi
  done < "$MAN"
  if [ "$claimed" -eq 0 ]; then
    run_duplicate_links
    echo "[NO_MORE_READY] $(date '+%F %T %Z')"
    break
  fi
done

echo "[END] $(date '+%F %T %Z') node=$NODE pid=$$"
