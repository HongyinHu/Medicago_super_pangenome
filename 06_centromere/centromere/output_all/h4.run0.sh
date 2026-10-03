#让背景更干净，可以把低于某个 log2 值的小波动设为 0。比如低于 0.5 的信号全部隐藏

SP=genome_474

BW=02_chipseq/genome_474/05_chipseq_signal/genome_474_CENH3_vs_Input.q0.rmdup.100000.smooth500000.log2.bw
CHROMSIZE=genome_474.chrom.sizes

BG=${BW}.bedGraph
BG_FILTER=${BW}.pc10.log2.min0.5.bedGraph

# BG=${SP}_CENH3_vs_Input.q0.rmdup.100k.smooth500k.pc10.log2.bedGraph
# BG_FILTER=${SP}_CENH3_vs_Input.q0.rmdup.100k.smooth500k.pc10.log2.min0.5.bedGraph


bigWigToBedGraph ${BW} ${BG}

awk 'BEGIN{OFS="\t"}
{
    if($4 < 0.5) $4 = 0;
    print
}' ${BG} > ${BG_FILTER}

bedGraphToBigWig \
  ${BG_FILTER} \
  ${CHROMSIZE} \
  ${BG_FILTER}.bw