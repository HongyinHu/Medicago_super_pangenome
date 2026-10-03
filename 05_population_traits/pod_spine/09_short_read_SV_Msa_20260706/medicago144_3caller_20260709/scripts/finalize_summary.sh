#!/usr/bin/env bash
set -euo pipefail
MED=path/to/project/N_4.pod_spiny/09_short_read_SV_Msa_20260706/medicago144_3caller_20260709
BCFTOOLS=path/to/home/anaconda3/envs/panpop/bin/bcftools
cd "$MED"
N=$(($(wc -l < input/medicago144_sample_bams.tsv)-1))
delly_done=$(find results/delly_per_sample -name '*.delly.all.done' | wc -l)
manta_done=$(find results/manta_per_sample -name '*.manta.done' | wc -l)
manta_fail=$(find results/manta_per_sample -name '*.manta.failed' | wc -l)
smoove_done=$(find results/smoove_per_sample -name '*.smoove.done' | wc -l)
smoove_fail=$(find results/smoove_per_sample -name '*.smoove.failed' | wc -l)
cons_done=$(find results/consensus_per_sample -name '*.consensus_2of3.done' | wc -l)
cons_vcf=$(find results/consensus_per_sample -name '*.shortread_3caller.consensus2.vcf.gz' | wc -l)
{
  echo -e "metric\tvalue"
  echo -e "expected_samples\t$N"
  echo -e "delly_done\t$delly_done"
  echo -e "manta_done\t$manta_done"
  echo -e "manta_failed_unavailable\t$manta_fail"
  echo -e "smoove_done\t$smoove_done"
  echo -e "smoove_failed_unavailable\t$smoove_fail"
  echo -e "consensus_done\t$cons_done"
  echo -e "consensus_vcf\t$cons_vcf"
  echo -e "time\t$(date '+%F %T %Z')"
} > summary/medicago144_3caller_progress.tsv
cp -f summary/medicago144_3caller_progress.tsv summary/progress.tsv
# per-sample consensus count
{
  echo -e "sample_id\tconsensus_records\tdelly_available\tmanta_available\tsmoove_available"
  while IFS=$'\t' read -r sid rest; do
    [ "$sid" = "sample_id" ] && continue
    donefile="results/consensus_per_sample/$sid/$sid.consensus_2of3.done"
    rec=NA; da=0; ma=0; sa=0
    [ -s "results/delly_per_sample/$sid/$sid.delly.all.vcf.gz" ] && da=1
    [ -s "results/manta_per_sample/$sid/results/variants/diploidSV.vcf.gz" ] && ma=1
    [ -s "results/smoove_per_sample/$sid/$sid-smoove.genotyped.vcf.gz" ] && sa=1
    if [ -s "$donefile" ]; then rec=$(awk -F'\t' '$1=="records"{print $2}' "$donefile"); [ -z "$rec" ] && rec=0; fi
    echo -e "$sid\t$rec\t$da\t$ma\t$sa"
  done < input/medicago144_sample_bams.tsv
} > summary/medicago144_consensus_record_counts.tsv
find logs -type f \( -name '*.err' -o -name '*.out' \) -size +0c -print0 | xargs -0 -r grep -IlE 'No space left|Disk quota exceeded|ENOSPC|OUT_OF_MEMORY|OOM|Killed|FATAL|ERROR|Exception|Traceback|not found|Permission' > summary/logs_with_errors.txt || true
cat summary/medicago144_3caller_progress.tsv
