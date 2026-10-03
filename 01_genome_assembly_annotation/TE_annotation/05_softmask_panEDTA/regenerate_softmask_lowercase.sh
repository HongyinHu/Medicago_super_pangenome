#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/N_1.EDTA_single
LAI=$BASE/07_LAI_panEDTA_IDfixed_fullpass
SOFT=$BASE/05_softmask_panEDTA
LOGDIR=$SOFT/logs

export PATH=path/to/home/anaconda3/envs/EDTA_env/bin:$PATH

GENOMES=(
  genome_395.fa genome_410.fa genome_436.fa genome_454.fa genome_457.fa
  genome_461.fa genome_468.fa genome_472.fa genome_474a.fa genome_474b.fa
  genome_482.fa genome_A17.fa genome_M22.fa genome_M46.fa genome_Mar.fa
  genome_Mpo.fa genome_Mru.fa genome_Msa1.fa genome_Msa2.fa genome_R108.fa
  genome_ZM4.fa
)

mkdir -p "$SOFT" "$LOGDIR" "$SOFT/bed"
exec > >(tee -a "$LOGDIR/regenerate_softmask_lowercase.log") 2>&1

ts() {
  date '+%F %T %Z'
}

out_to_bed() {
  local out=$1
  local bed=$2
  awk 'BEGIN{OFS="\t"}
    $1 ~ /^[0-9]+$/ && $6 ~ /^[0-9]+$/ && $7 ~ /^[0-9]+$/ {
      start=$6-1
      end=$7
      if (start < 0) start=0
      if (end < start) {
        tmp=start
        start=end-1
        end=tmp+1
      }
      print $5, start, end
    }' "$out" | sort -k1,1 -k2,2n | bedtools merge -i - > "$bed"
}

base_counts() {
  local fa=$1
  awk '
    /^>/ {next}
    {
      s=$0
      lower += gsub(/[acgt]/, "", s)
      s=$0
      upper += gsub(/[ACGT]/, "", s)
      s=$0
      n += gsub(/[Nn]/, "", s)
    }
    END {printf "%d\t%d\t%d\n", upper+0, lower+0, n+0}
  ' "$fa"
}

echo "START $(ts) host=$(hostname)"
command -v bedtools

old="$SOFT/fasta"
if [ -d "$old" ]; then
  stamp=$(date +%Y%m%d_%H%M%S)
  mv "$old" "$SOFT/fasta_partial_or_previous_backup_$stamp"
  echo "BACKUP $(ts) moved $old to $SOFT/fasta_partial_or_previous_backup_$stamp"
fi
mkdir -p "$SOFT/fasta"

printf 'sample\toriginal_genome\trmout_idfixed\tbed\tsoftmask_fasta\tfasta_bytes\tupper_acgt\tlower_acgt\tn_count\tstatus\n' > "$SOFT/manifest.tsv.tmp"

failed=0
for g in "${GENOMES[@]}"; do
  sample=${g%.fa}
  mod=$LAI/$sample/$g.mod
  rmout=$LAI/$sample/$g.mod.panEDTA.IDfixed.out
  bed=$SOFT/bed/$sample.panEDTA.repeatmasker.bed
  fa=$SOFT/fasta/$sample.panEDTA.softmasked.fa
  status=ok

  if [ ! -s "$mod" ] || [ ! -s "$rmout" ]; then
    echo "ERROR $(ts) $sample missing input mod=$mod rmout=$rmout"
    failed=$((failed + 1))
    status=missing_input
    printf '%s\t%s\t%s\t%s\t%s\t0\t0\t0\t0\t%s\n' "$sample" "$mod" "$rmout" "$bed" "$fa" "$status" >> "$SOFT/manifest.tsv.tmp"
    continue
  fi

  echo "MASK $(ts) $sample"
  out_to_bed "$rmout" "$bed.tmp"
  mv "$bed.tmp" "$bed"
  awk '/^>/ {print; next} {print toupper($0)}' "$mod" > "$fa.input_upper.tmp"
  bedtools maskfasta -soft -fi "$fa.input_upper.tmp" -bed "$bed" -fo "$fa.tmp"
  mv "$fa.tmp" "$fa"
  rm -f "$fa.input_upper.tmp"

  read -r upper lower ncount < <(base_counts "$fa")
  bytes=$(stat -c '%s' "$fa")
  if [ "$lower" -le 0 ]; then
    status=no_lowercase
    failed=$((failed + 1))
    echo "ERROR $(ts) $sample generated no lowercase bases"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$mod" "$rmout" "$bed" "$fa" "$bytes" "$upper" "$lower" "$ncount" "$status" >> "$SOFT/manifest.tsv.tmp"
done

mv "$SOFT/manifest.tsv.tmp" "$SOFT/manifest.tsv"
find "$SOFT/fasta" -maxdepth 1 -type f -name '*.panEDTA.softmasked.fa' | wc -l > "$SOFT/softmask.count"
echo "DONE $(ts) count=$(cat "$SOFT/softmask.count") failed=$failed manifest=$SOFT/manifest.tsv"
exit "$failed"
