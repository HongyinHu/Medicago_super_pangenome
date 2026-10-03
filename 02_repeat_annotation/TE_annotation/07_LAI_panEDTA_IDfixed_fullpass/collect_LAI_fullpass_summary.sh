#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/N_1.EDTA_single
OUTDIR=$BASE/07_LAI_panEDTA_IDfixed_fullpass
LOGROOT=$OUTDIR/logs
SUMMARY=$OUTDIR/LAI_summary_fullpass.tsv
QUEUE=$OUTDIR/genome.queue

printf 'genome\tLAI_file\tLAI_value\tintact_pass_list\tpanEDTA_IDfixed_out\tstatus\tnote\n' > "$SUMMARY.tmp"
while IFS= read -r genome; do
  [ -n "$genome" ] || continue
  sample=${genome%.fa}
  row="$LOGROOT/status/$sample.tsv"
  lai="$OUTDIR/$sample/$genome.mod.panEDTA.IDfixed.out.LAI"
  if [ -s "$lai" ]; then
    value=$(awk 'BEGIN{v="NA"} $1=="whole_genome"{v=$7} END{print v}' "$lai" 2>/dev/null || echo NA)
    intact=$(readlink -f "$OUTDIR/$sample/$genome.mod.pass.list" 2>/dev/null || echo NA)
    panout=$(readlink -f "$OUTDIR/$sample/$genome.mod.panEDTA.IDfixed.out" 2>/dev/null || echo NA)
    printf '%s\t%s\t%s\t%s\t%s\tdone\tfinal_present\n' "$sample" "$lai" "$value" "$intact" "$panout" >> "$SUMMARY.tmp"
  elif [ -s "$row" ]; then
    cat "$row" >> "$SUMMARY.tmp"
  else
    printf '%s\t%s\tNA\tNA\tNA\tmissing\tno_status\n' "$sample" "$lai" >> "$SUMMARY.tmp"
  fi
done < "$QUEUE"
mv "$SUMMARY.tmp" "$SUMMARY"
echo "[$(date)] wrote $SUMMARY" | tee -a "$LOGROOT/collect.log"
