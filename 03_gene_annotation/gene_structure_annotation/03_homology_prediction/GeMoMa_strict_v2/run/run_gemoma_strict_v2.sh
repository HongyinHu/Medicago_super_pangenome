#!/usr/bin/env bash
set -euo pipefail
SAMPLE=${1:?sample}
THREADS=${THREADS:-8}
HEAP=${HEAP:-300g}
BASE=path/to/project/N_1.coding_gene_anno
GM=$BASE/03_homology_prediction/GeMoMa_strict_v2
REFS=$GM/refs/reference_inputs/related_gemoma_refs.tsv
MANIFEST=$BASE/01_genome_versions/genome_manifest.tsv
LOCKROOT=$GM/locks
LOGROOT=$GM/logs
ASCII_REF_DIR=$GM/refs/ascii_for_direct_java
JAVA=path/to/home/anaconda3/envs/software/lib/jvm/bin/java
JAR=path/to/home/anaconda3/envs/software/share/gemoma-1.9-0/GeMoMa-1.9.jar
BIN=path/to/home/anaconda3/envs/software/bin
mkdir -p "$LOCKROOT" "$LOGROOT" "$ASCII_REF_DIR" "$GM/failed_attempts"
export LANG=${LANG:-C.UTF-8}
export LC_ALL=${LC_ALL:-C.UTF-8}
host=$(hostname)
stamp=$(date +%Y%m%d_%H%M%S)
lock=$LOCKROOT/${SAMPLE}.lock
outdir=$GM/$SAMPLE
done_file=$outdir/gemoma_strict_v2.done
fail_file=$outdir/gemoma_strict_v2.failed
log=$LOGROOT/strict_v2_${SAMPLE}_${host}_${stamp}.log
exec > >(tee -a "$log") 2>&1

echo "[$(date)] START GeMoMa strict_v2 sample=$SAMPLE host=$host threads=$THREADS heap=$HEAP"
if [ -s "$done_file" ]; then echo "[$(date)] SKIP already done $done_file"; exit 0; fi
if ! mkdir "$lock" 2>/dev/null; then echo "[$(date)] SKIP lock exists $lock"; exit 9; fi
echo "host=$host" > "$lock/owner"
echo "pid=$$" >> "$lock/owner"
echo "started=$(date '+%F %T %Z')" >> "$lock/owner"
echo "profile=strict_v2 self_reference_excluded GeMoMa.p=1 GeMoMa.e=1E-5 GAF.k=YES evidence>1" >> "$lock/owner"
cleanup(){ rc=$?; if [ "$rc" -eq 0 ]; then rm -rf "$lock"; else echo "failed_rc=$rc $(date '+%F %T %Z')" >> "$lock/owner" 2>/dev/null || true; fi; exit "$rc"; }
trap cleanup EXIT

target=$(awk -v s="$SAMPLE" 'BEGIN{FS="\t"} NR>1 && $1==s{print $3; exit}' "$MANIFEST")
if [ -z "$target" ] || [ ! -s "$target" ]; then
  mkdir -p "$outdir"
  echo -e "failed\tmissing_target\t$(date)\t$host" > "$fail_file"
  echo "[$(date)] ERROR missing target genome: $target"
  exit 2
fi
if [ -d "$outdir" ] && [ ! -s "$done_file" ]; then
  archive=$GM/failed_attempts/${SAMPLE}_strict_v2_archive_${stamp}
  mv "$outdir" "$archive"
  echo "[$(date)] archived old strict outdir -> $archive"
fi
mkdir -p "$outdir"
refs=()
while IFS=$'\t' read -r ref genome gff weight status; do
  [ "$ref" = "ref_id" ] && continue
  [ -n "$ref" ] || continue
  case "$SAMPLE:$ref" in
    genome_A17:A17*|genome_R108:R108*|genome_Mar:Mar|genome_Mpo:Mpo|genome_Mru:Mru|genome_Msa1:Msa|genome_Msa2:Msa) continue ;;
  esac
  safe_ref=$(printf '%s' "$ref" | tr -c 'A-Za-z0-9_.-' '_')
  ascii_genome=$ASCII_REF_DIR/${safe_ref}.genome.fa
  ascii_gff=$ASCII_REF_DIR/${safe_ref}.annotation.gff
  ln -sfn "$genome" "$ascii_genome"
  ln -sfn "$gff" "$ascii_gff"
  refs+=(s=own i="$ref" a="$ascii_gff" g="$ascii_genome" w="$weight")
done < "$REFS"

echo "[$(date)] RUN strict_v2 sample=$SAMPLE nref=$(( ${#refs[@]} / 5 )) target=$target"
set +e
"$JAVA" -Xms4g -Xmx"$HEAP" -jar "$JAR" CLI GeMoMaPipeline \
  t="$target" "${refs[@]}" \
  tblastn=false m="$BIN" b="$BIN" sc=false debug=false p=true pc=true o=true \
  GeMoMa.e=1E-5 GeMoMa.p=1 GAF.k=YES \
  "GAF.f=start=='M' and stop=='*' and score/aa>=0.75 and evidence>1" \
  "GAF.a=sumWeight>1" GAF.mnotpg=2 \
  AnnotationFinalizer.r=SIMPLE AnnotationFinalizer.p=${SAMPLE}_SG \
  outdir="$outdir" threads="$THREADS" \
  > "$outdir/gemoma.stdout.log" 2> "$outdir/gemoma.stderr.log"
rc=$?
set -e
if [ "$rc" -eq 0 ] && [ -s "$outdir/final_annotation.gff" ]; then
  genes=$(awk '$3=="gene"{n++} END{print n+0}' "$outdir/final_annotation.gff")
  ev1=$(awk 'BEGIN{genes=0; ev1=0} $3=="gene"{genes++; if($0 ~ /maxEvidence=1(;|$)/) ev1++} END{if(genes>0) printf "%.1f", ev1*100/genes; else print "NA"}' "$outdir/final_annotation.gff")
  echo -e "ok\t$(date)\t$host\t$THREADS\theap=$HEAP\tgenes=$genes\tmaxEvidence1_pct=$ev1\tprofile=strict_v2\tself_reference_excluded" > "$done_file"
  rm -f "$fail_file"
  echo "[$(date)] DONE strict_v2 $SAMPLE genes=$genes maxEvidence1_pct=$ev1"
else
  echo -e "failed\trc=$rc\t$(date)\t$host\theap=$HEAP\tprofile=strict_v2" > "$fail_file"
  echo "[$(date)] FAILED strict_v2 $SAMPLE rc=$rc"
  tail -80 "$outdir/gemoma.stderr.log" || true
  tail -80 "$outdir/gemoma.stdout.log" || true
  exit "$rc"
fi
