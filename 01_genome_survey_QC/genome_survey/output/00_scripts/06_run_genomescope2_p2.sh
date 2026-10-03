#!/usr/bin/env bash
set -euo pipefail
OUT=path/to/project/N_2.DNA_survy/output
R=path/to/home/anaconda3/envs/R/bin/Rscript
GS="$OUT/software/genomescope2.0/genomescope.R"
RES="$OUT/04_genomescope2_p2"
mkdir -p "$RES"
STATUS="$RES/genomescope2_p2_status.tsv"
printf "time\tsample\tstatus\tnote\n" > "$STATUS"
for hist in "$OUT"/02_jellyfish/k21/*.histo; do
  sample=$(basename "$hist" .k21.histo)
  od="$RES/$sample"
  mkdir -p "$od"
  printf "%s\t%s\tstart\t%s\n" "$(date '+%F %T')" "$sample" "$hist" >> "$STATUS"
  if "$R" "$GS" -i "$hist" -o "$od" -p 2 -k 21 -n "$sample" --json_report > "$od/run.log" 2>&1; then
    for png in "$od"/*.png; do
      [ -f "$png" ] || continue
      convert "$png" "${png%.png}.pdf" >/dev/null 2>&1 || true
    done
    printf "%s\t%s\tdone\t%s\n" "$(date '+%F %T')" "$sample" "$od" >> "$STATUS"
  else
    printf "%s\t%s\tfailed\t%s/run.log\n" "$(date '+%F %T')" "$sample" "$od" >> "$STATUS"
  fi
done
python3 - <<'PY'
from pathlib import Path
import re
out = Path('path/to/project/N_2.DNA_survy/output/04_genomescope2_p2')
rows=[]
for summ in sorted(out.glob('*/**/*_summary.txt')):
    sample=summ.parent.name
    vals={'sample':sample,'summary':str(summ)}
    text=summ.read_text(errors='replace')
    for line in text.splitlines():
        if line.startswith('Genome Haploid Length'):
            nums=re.findall(r'[0-9,]+ bp', line)
            vals['haploid_min_bp']=nums[0].replace(',','').replace(' bp','') if len(nums)>0 else ''
            vals['haploid_max_bp']=nums[1].replace(',','').replace(' bp','') if len(nums)>1 else ''
        elif line.startswith('Genome Repeat Length'):
            nums=re.findall(r'[0-9,]+ bp', line)
            vals['repeat_min_bp']=nums[0].replace(',','').replace(' bp','') if len(nums)>0 else ''
            vals['repeat_max_bp']=nums[1].replace(',','').replace(' bp','') if len(nums)>1 else ''
        elif line.startswith('Genome Unique Length'):
            nums=re.findall(r'[0-9,]+ bp', line)
            vals['unique_min_bp']=nums[0].replace(',','').replace(' bp','') if len(nums)>0 else ''
            vals['unique_max_bp']=nums[1].replace(',','').replace(' bp','') if len(nums)>1 else ''
        elif line.startswith('Model Fit'):
            nums=re.findall(r'[0-9.]+%', line)
            vals['model_fit_min_pct']=nums[0].replace('%','') if len(nums)>0 else ''
            vals['model_fit_max_pct']=nums[1].replace('%','') if len(nums)>1 else ''
        elif line.startswith('Read Error Rate'):
            nums=re.findall(r'[0-9.]+%', line)
            vals['read_error_pct']=nums[0].replace('%','') if nums else ''
        elif line.startswith('Heterozygous'):
            nums=re.findall(r'[0-9.]+%', line)
            vals['heterozygous_min_pct']=nums[0].replace('%','') if len(nums)>0 else ''
            vals['heterozygous_max_pct']=nums[1].replace('%','') if len(nums)>1 else ''
    rows.append(vals)
fields=['sample','haploid_min_bp','haploid_max_bp','repeat_min_bp','repeat_max_bp','unique_min_bp','unique_max_bp','heterozygous_min_pct','heterozygous_max_pct','model_fit_min_pct','model_fit_max_pct','read_error_pct','summary']
with open(out/'genomescope2_p2_summary.tsv','w') as fh:
    fh.write('\t'.join(fields)+'\n')
    for r in rows:
        fh.write('\t'.join(r.get(f,'') for f in fields)+'\n')
PY
