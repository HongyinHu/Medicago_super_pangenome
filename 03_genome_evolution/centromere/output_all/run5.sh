SP=genome_474

bedtools makewindows \
  -g ${SP}.chrom.sizes \
  -w 100000 \
  > 05_centromere_define/${SP}.100kb.window.bed

multiBigwigSummary BED-file \
  --BED 05_centromere_define/${SP}.100kb.window.bed \
  -b 04_cenh3_signal/${SP}.CENH3_vs_Input.q20.10k.log2.bw \
  -o 05_centromere_define/${SP}.CENH3.100kb.npz \
  --outRawCounts 05_centromere_define/${SP}.CENH3.100kb.signal.tab
