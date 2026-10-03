#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/N_1.coding_gene_anno
PASA_IN=$BASE/04_transcript_prediction/PASA_inputs
PB=$PASA_IN/pacbio_bam_processing
OUT=$PASA_IN/combined
MAN=$PASA_IN/manifests/pasa_transcripts_manifest.tsv

mkdir -p "$OUT" "$PASA_IN/manifests"
printf 'sample\tcombined_fasta\tn_total\tn_illumina\tn_isoseq_pacbio\tn_isoseq_existing\tstatus\tsources\n' > "$MAN"

clean_one() {
  local sample="$1" label="$2" infile="$3"
  awk -v sample="$sample" -v label="$label" '
    BEGIN{RS=">"; ORS=""; n=0}
    NR>1{
      split($0,a,"\n");
      h=a[1];
      gsub(/[^A-Za-z0-9_.:|=-]/,"_",h);
      prefix=sample "|" label "|";
      if(index(h,prefix)==1){
        rest=substr(h,length(prefix)+1);
        split(rest,b,"|");
        if(b[1] ~ /^[0-9]+$/){
          h=substr(rest,length(b[1])+2);
        }
      }
      if(length(h)==0){h="seq"}
      seq="";
      for(i=2;i<=length(a);i++) seq=seq a[i];
      gsub(/[[:space:]]/,"",seq);
      gsub(/[^A-Za-z]/,"",seq);
      seq=toupper(seq);
      if(length(seq)>=100){
        n++;
        print ">"sample"|"label"|"n"|"h"\n";
        for(i=1;i<=length(seq);i+=80) print substr(seq,i,80)"\n";
      }
    }' "$infile"
}

tail -n +2 "$BASE/01_genome_versions/genome_manifest.tsv" | while IFS=$'\t' read -r sample soft unmask rest; do
  [ -n "$sample" ] || continue
  sources=()
  labels=()

  illum="$PASA_IN/stringtie/$sample/${sample}.illumina_stringtie.transcripts.fa"
  pb_iso="$PB/final/${sample}.isoseq_transcripts.fa"
  existing_iso="$PASA_IN/isoseq_existing_fasta/$sample/${sample}.isoseq_existing.transcripts.fa"

  [ -s "$illum" ] && { sources+=("$illum"); labels+=("illumina_stringtie"); }
  if [ -s "$pb_iso" ]; then
    sources+=("$pb_iso"); labels+=("pacbio_isoseq")
  elif [ -s "$existing_iso" ]; then
    sources+=("$existing_iso"); labels+=("isoseq_existing")
  fi

  out="$OUT/${sample}.pasa_transcripts.fa"
  tmp="$out.tmp"
  : > "$tmp"

  if [ "${#sources[@]}" -eq 0 ]; then
    rm -f "$tmp"
    printf '%s\t%s\t0\t0\t0\t0\tmissing\t\n' "$sample" "$out" >> "$MAN"
    continue
  fi

  nillum=0
  npb=0
  nexisting=0
  for i in "${!sources[@]}"; do
    src="${sources[$i]}"
    label="${labels[$i]}"
    clean_one "$sample" "$label" "$src" >> "$tmp"
    count=$(grep -c '^>' "$src" 2>/dev/null || echo 0)
    case "$label" in
      illumina_stringtie) nillum=$count ;;
      pacbio_isoseq) npb=$count ;;
      isoseq_existing) nexisting=$count ;;
    esac
  done

  awk 'BEGIN{RS=">"; ORS=""}
    NR>1{
      split($0,a,"\n");
      h=a[1];
      seq="";
      for(i=2;i<=length(a);i++) seq=seq a[i];
      gsub(/[[:space:]]/,"",seq);
      if(length(seq)>=100 && !seen[seq]++){
        print ">"h"\n";
        for(i=1;i<=length(seq);i+=80) print substr(seq,i,80)"\n";
      }
    }' "$tmp" > "$out"
  rm -f "$tmp"
  rm -f "$out.fai"

  ntot=$(grep -c '^>' "$out" 2>/dev/null || echo 0)
  if [ "$ntot" -gt 0 ]; then
    status=ok
  else
    status=empty
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$out" "$ntot" "$nillum" "$npb" "$nexisting" "$status" "$(IFS=,; echo "${sources[*]}")" >> "$MAN"
done
