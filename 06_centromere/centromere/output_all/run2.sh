SP=genome_474

bamCompare \
  -b1 04_cenh3_signal/${SP}.CENH3.q20.bam \
  -b2 04_cenh3_signal/${SP}.Input.q20.bam \
  --operation log2 \
  --binSize 5000 \
  --normalizeUsing CPM \
  --scaleFactorsMethod None \
  -o 04_cenh3_signal/${SP}.CENH3_vs_Input.q20.10k.log2.bw

bamCompare \
  -b1 04_cenh3_signal/${SP}.CENH3.relaxed.bam \
  -b2 04_cenh3_signal/${SP}.Input.relaxed.bam \
  --operation log2 \
  --binSize 10000 \
  --normalizeUsing CPM \
  --scaleFactorsMethod None \
  -o 04_cenh3_signal/${SP}.CENH3_vs_Input.relaxed.10k.log2.bw
