#!/usr/bin/env bash
set -euo pipefail

BASE=path/to/project/N_1.coding_gene_anno
REF=$BASE/00_inputs/related_species_reference
AB=$BASE/02_abinitio_prediction
OUTDIR=$AB/inputs
OUT=$OUTDIR/related_species.clean.pep.fa
MAN=$OUTDIR/related_species.clean.pep.manifest.tsv

mkdir -p "$OUTDIR"
printf 'source\tseq_count\tstatus\n' > "$MAN"
tmp="$OUT.tmp"
: > "$tmp"

find "$REF" -type f \( -iname '*.pep' -o -iname '*.pep.fa' -o -iname '*.faa' -o -iname '*.protein*.fa' -o -path '*/derived_pep/*.fa' \) | sort | while read -r f; do
  n=$(grep -c '^>' "$f" 2>/dev/null || echo 0)
  if [ "$n" -gt 0 ]; then
    printf '%s\t%s\tok\n' "$f" "$n" >> "$MAN"
    awk -v src="$(basename "$f")" '
      BEGIN{RS=">"; ORS=""}
      NR>1{
        split($0,a,"\n");
        h=a[1];
        gsub(/[^A-Za-z0-9_.:|=-]/,"_",h);
        seq="";
        for(i=2;i<=length(a);i++) seq=seq a[i];
        gsub(/[[:space:]]/,"",seq);
        seq=toupper(seq);
        gsub(/[^A-Z*]/,"X",seq);
        gsub(/\*/,"",seq);
        if(length(seq)>=50){
          print ">"src"|"h"\n";
          for(i=1;i<=length(seq);i+=80) print substr(seq,i,80)"\n";
        }
      }' "$f" >> "$tmp"
  fi
done

awk 'BEGIN{RS=">"; ORS=""}
  NR>1{
    split($0,a,"\n");
    h=a[1];
    seq="";
    for(i=2;i<=length(a);i++) seq=seq a[i];
    gsub(/[[:space:]]/,"",seq);
    if(length(seq)>=50 && !seen[seq]++){
      print ">"h"\n";
      for(i=1;i<=length(seq);i+=80) print substr(seq,i,80)"\n";
    }
  }' "$tmp" > "$OUT"
rm -f "$tmp"
printf 'TOTAL\t%s\t%s\n' "$(grep -c '^>' "$OUT" 2>/dev/null || echo 0)" "$OUT" >> "$MAN"
