SP=genome_474

GSIZE=$(awk '{s+=$2}END{print s}' ${SP}.chrom.sizes)

path/to/home/anaconda3/envs/CENH3_env/bin/macs3 callpeak \
  -t 04_cenh3_signal/${SP}.CENH3.q20.bam \
  -c 04_cenh3_signal/${SP}.Input.q20.bam \
  -f BAMPE \
  -g $GSIZE \
  --broad \
  --broad-cutoff 0.1 \
  --keep-dup all \
  -n ${SP}.CENH3.q20 \
  --outdir 04_cenh3_signal/${SP}.macs2

echo "1. finshed macs3"

#CENH3_DOMAIN=04_cenh3_signal/${SP}.macs2/${SP}.CENH3.q20_peaks.broadPeak.bed

cut -f1-3 \
  04_cenh3_signal/${SP}.macs2/${SP}.CENH3.q20_peaks.broadPeak \
  > 04_cenh3_signal/${SP}.macs2/${SP}.CENH3.q20_peaks.broadPeak.bed

# 在脚本开头或者该命令前添加环境变量
export LC_ALL=C

sort -k1,1 -k2,2n 04_cenh3_signal/${SP}.macs2/${SP}.CENH3.q20_peaks.broadPeak.bed | \
   bedtools merge -i -  -d 100000 > 04_cenh3_signal/${SP}.macs2/${SP}.CENH3.q20_peaks.broadPeak.merge100k.bed

echo "finshed merge"


