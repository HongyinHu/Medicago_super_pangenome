#!/usr/bin/env bash
set -euo pipefail
BASE=${BASE:-path/to/project/N_3.call_SV}
OUT=${OUT:-$BASE/02.assembly_svgap_dualref}
SVGAP=${SVGAP:-$OUT/software/SVGAP}
REF_RUN=${1:?usage: run_svgap_ref_stages.sh Msa_ref|R108_ref [THREADS]}
THREADS=${2:-12}
case "$REF_RUN" in
  Msa_ref) REF=Msa ;;
  R108_ref) REF=R108 ;;
  *) echo "bad REF_RUN=$REF_RUN" >&2; exit 2 ;;
esac
TASKS=$OUT/01_metadata/svgap_wga_tasks.tsv
need=$(awk -F '\t' -v r="$REF_RUN" 'NR>1 && $2==r{n++} END{print n+0}' "$TASKS")
done_n=$(awk -F '\t' -v r="$REF_RUN" -v out="$OUT" 'NR>1 && $2==r{f=out"/status/wga/"$1".done"; if((getline x < f)>=0){n++}; close(f)} END{print n+0}' "$TASKS")
if [ "$done_n" -lt "$need" ] && [ "${ALLOW_PARTIAL:-0}" != "1" ]; then
  echo "WGA_NOT_COMPLETE ref_run=$REF_RUN done=$done_n need=$need" >&2
  exit 3
fi
cd "$OUT/$REF_RUN"
mkdir -p chainnet
log="$OUT/logs/${REF_RUN}.svgap_stages.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$log") 2>&1
echo "START_SVGAP_STAGES ref_run=$REF_RUN ref=$REF threads=$THREADS time=$(date '+%F %T %Z')"
perl "$SVGAP/1_Convert2Axt.pl" --ali minimap2 --input "$OUT/$REF_RUN/alignment" --wk "$OUT/$REF_RUN/chainnet" --tname "$REF" --t "$THREADS"
perl "$SVGAP/2_ChainNetSyn.pl" --gd "$OUT/genome" --ad "$OUT/$REF_RUN/chainnet/Target_$REF" --lst "$OUT/genome/g.lst" --sing --syn synnet
perl "$SVGAP/3_SynNetFilter.pl" --synnet synnet --gd "$OUT/genome" --chain "$OUT/$REF_RUN/chainnet/Target_$REF"
perl "$SVGAP/4_PairwiseSV.pl" --syndir CleanSynNet --gd "$OUT/genome" --outdir CleanSV --t "$THREADS"
perl "$SVGAP/5_Combined.pl" --sv CleanSV --refname "$REF" --outdir CombinedSV
touch "$OUT/status/${REF_RUN}.svgap_stages.done"
echo "DONE_SVGAP_STAGES ref_run=$REF_RUN time=$(date '+%F %T %Z')"
