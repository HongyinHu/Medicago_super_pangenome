#!/usr/bin/env bash
set -euo pipefail
sample=${1:?sample}
BASE=path/to/project/N_1.coding_gene_anno
E="$BASE/05_EVM_integration"
row=$(awk -F'\t' -v s="$sample" 'NR>1 && $1==s && $7=="ready"{print}' "$E/ready_samples.tsv")
if [ -z "$row" ]; then echo "sample_not_ready_or_missing_in_manifest: $sample" >&2; exit 3; fi
IFS=$'\t' read -r sample genome braker_gff braker_src gemoma_gff pasa_gff status expected_signature <<< "$row"
[[ "$gemoma_gff" == */GeMoMa_strict_v2/* ]] || { echo "non_strict_v2_gemoma_input: $gemoma_gff" >&2; exit 5; }

current_signature=$(
  for path in "$genome" "$braker_gff" "$gemoma_gff" "$pasa_gff"; do
    [ -s "$path" ] || { echo "missing_input_for_signature: $path" >&2; exit 6; }
    stat -Lc '%n\t%s\t%Y' "$path"
  done | sha256sum | awk '{print $1}'
)
[ "$current_signature" = "$expected_signature" ] || { echo "stale_ready_manifest expected=$expected_signature current=$current_signature" >&2; exit 7; }

out="$E/work/$sample"
data="$out/00_inputs"
tmp="$out/00_inputs.tmp.$$"
rm -rf "$tmp"
mkdir -p "$tmp"
data="$tmp"
ln -sfn "$genome" "$data/genome.fa"
python3 - "$braker_gff" "$data/braker_predictions.gff3" BRAKER3 <<'PY'
import sys
inp,out,src=sys.argv[1:4]
keep={"gene","mRNA","transcript","exon","CDS","start_codon","stop_codon"}
with open(inp) as f, open(out,'w') as w:
    w.write('##gff-version 3\n')
    for line in f:
        if line.startswith('#'): continue
        parts=line.rstrip('\n').split('\t')
        if len(parts)<9 or parts[2] not in keep: continue
        if parts[2]=='transcript': parts[2]='mRNA'
        parts[1]=src
        w.write('\t'.join(parts)+'\n')
PY
python3 - "$gemoma_gff" "$data/gemoma_predictions.gff3" GeMoMa <<'PY'
import sys
inp,out,src=sys.argv[1:4]
keep={"gene","mRNA","transcript","exon","CDS","start_codon","stop_codon"}
with open(inp) as f, open(out,'w') as w:
    w.write('##gff-version 3\n')
    for line in f:
        if line.startswith('#'): continue
        parts=line.rstrip('\n').split('\t')
        if len(parts)<9 or parts[2] not in keep: continue
        if parts[2]=='transcript': parts[2]='mRNA'
        parts[1]=src
        w.write('\t'.join(parts)+'\n')
PY
{ grep '^#' "$data/braker_predictions.gff3" | head -1; grep -v '^#' "$data/braker_predictions.gff3"; grep -v '^#' "$data/gemoma_predictions.gff3"; } > "$data/gene_predictions.gff3"
python3 - "$pasa_gff" "$data/transcript_alignments.gff3" PASA <<'PY'
import sys
inp,out,src=sys.argv[1:4]
keep={"gene","mRNA","transcript","exon","CDS","cDNA_match","match","match_part"}
with open(inp) as f, open(out,'w') as w:
    w.write('##gff-version 3\n')
    for line in f:
        if line.startswith('#'): continue
        parts=line.rstrip('\n').split('\t')
        if len(parts)<9 or parts[2] not in keep: continue
        if parts[2]=='transcript': parts[2]='mRNA'
        parts[1]=src
        w.write('\t'.join(parts)+'\n')
PY
if [ "$braker_src" = "BRAKER3_RNA_protein" ]; then
  braker_weight=5
  gemoma_weight=5
else
  braker_weight=3
  gemoma_weight=6
fi
printf 'ABINITIO_PREDICTION\tBRAKER3\t%s\nOTHER_PREDICTION\tGeMoMa\t%s\nTRANSCRIPT\tPASA\t10\n' \
  "$braker_weight" "$gemoma_weight" > "$data/weights.txt"
for f in "$data/genome.fa" "$data/gene_predictions.gff3" "$data/transcript_alignments.gff3" "$data/weights.txt"; do
  [ -s "$f" ] || { echo "missing_or_empty $f" >&2; exit 4; }
done
braker_genes=$(awk -F'\t' '$0!~/^#/ && $2=="BRAKER3" && $3=="gene"{n++} END{print n+0}' "$data/gene_predictions.gff3")
gemoma_genes=$(awk -F'\t' '$0!~/^#/ && $2=="GeMoMa" && $3=="gene"{n++} END{print n+0}' "$data/gene_predictions.gff3")
pasa_matches=$(awk -F'\t' '$0!~/^#/ && $2=="PASA" && ($3=="cDNA_match" || $3=="match"){n++} END{print n+0}' "$data/transcript_alignments.gff3")
[ "$braker_genes" -gt 0 ] && [ "$gemoma_genes" -gt 0 ] && [ "$pasa_matches" -gt 0 ] || {
  echo "invalid_evidence_counts braker=$braker_genes gemoma=$gemoma_genes pasa=$pasa_matches" >&2
  exit 8
}
cat > "$data/provenance.tsv" <<EOF
field	value
sample	$sample
input_signature	$current_signature
genome	$genome
braker_source	$braker_src
braker_gff3	$braker_gff
braker_genes	$braker_genes
gemoma_gff	$gemoma_gff
gemoma_genes	$gemoma_genes
pasa_gff3	$pasa_gff
pasa_matches	$pasa_matches
braker_weight	$braker_weight
gemoma_weight	$gemoma_weight
pasa_weight	10
prepared	$(date '+%F %T %Z')
EOF
rm -rf "$out/00_inputs"
mv "$tmp" "$out/00_inputs"
data="$out/00_inputs"
mkdir -p "$E/status"
echo "prepared $(date '+%F %T %Z') braker_source=$braker_src signature=$current_signature braker_genes=$braker_genes gemoma_genes=$gemoma_genes pasa_matches=$pasa_matches weights=BRAKER3:$braker_weight,GeMoMa:$gemoma_weight,PASA:10" > "$E/status/$sample.prepared"
