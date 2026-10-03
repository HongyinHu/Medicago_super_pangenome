#!/usr/bin/env bash
# Intra-genome Ks distribution of M. fischeriana with the NEW assembly (genome_M46_2),
# same pipeline and parameters as 5.WGD_2/output/X_genome_M46_1 (old assembly).
set -euo pipefail
W=path/to/project/5.WGD_2/output
OLD=$W/X_genome_M46_1
NEW=$W/X_genome_M46_2
SRC=path/to/project/37.karyotype_reconstruction/output_Kary/genome_M46_2
# the original CEGMA-bundled blast 2.10.0 no longer exists; biosofeware has the same
# blast 2.10.0+ and WGDI 0.5.7 (the version recorded for the earlier analyses)
BL=path/to/home/anaconda3/envs/biosofeware/bin
WGDI=path/to/home/anaconda3/envs/biosofeware/bin/wgdi
[ -x $WGDI ] || WGDI=path/to/home/anaconda3/envs/wgdi/bin/wgdi
echo "blastp: $($BL/blastp -version | head -1); wgdi: $WGDI"
mkdir -p $NEW && cd $NEW
rm -f status.done status.failed
trap 'rc=$?; [ $rc -ne 0 ] && touch status.failed; exit $rc' EXIT
{
echo "START $(date '+%F %T') $(hostname)"
cp $SRC/genome_M46_2.gff Mfi.gff
cp $SRC/genome_M46_2.lens Mfi.lens
cp $SRC/genome_M46_2.cds Mfi.cds
cp $SRC/genome_M46_2.pep Mfi.pep
echo "genes: $(wc -l < Mfi.gff) chromosomes: $(wc -l < Mfi.lens)"
# total.conf: identical to the old run, only the directory changes
sed "s#path/to/projects_all/pan_genome/pan_genome_construction/5.WGD_2/output/X_genome_M46_1#$NEW#g" \
  $OLD/total.conf > total.conf
grep -c X_genome_M46_2 total.conf
$BL/makeblastdb -in Mfi.pep -dbtype prot > /dev/null
$BL/blastp -num_threads 64 -db Mfi.pep -query Mfi.pep -outfmt 6 -evalue 1e-5 -num_alignments 20 -out Mfi.blastp.txt
echo "BLAST_DONE $(date '+%F %T')"
for step in -d -icl -ks -bi -c -bk -kp -pf; do
  echo "WGDI $step $(date '+%F %T')"; $WGDI $step total.conf
done
echo "END $(date '+%F %T')"
} > run_M46_2.log 2>&1
touch status.done
