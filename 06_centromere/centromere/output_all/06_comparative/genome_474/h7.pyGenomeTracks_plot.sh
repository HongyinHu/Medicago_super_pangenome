SP=genome_474

CHROMSIZE=../../${SP}.chrom.sizes

while read chr len
do
    echo "Plot ${chr}"

    pyGenomeTracks \
      --tracks tracks.ini \
      --region ${chr}:1-${len} \
      --outFileName ${SP}.${chr}_track.pdf \
      --width 32 \
      --dpi 300

done < ${CHROMSIZE}