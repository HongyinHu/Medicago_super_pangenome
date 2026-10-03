#环境CENH3_env

SP=genome_474

GENOME=00_genome/${SP}.fa
THREADS=32

samtools faidx $GENOME
cut -f1,2 ${GENOME}.fai > ${SP}.chrom.sizes

bwa-mem2 index $GENOME

bwa-mem2 mem -t $THREADS $GENOME \
	  02_chipseq/${SP}_CENH3_R1.fq.gz 02_chipseq/${SP}_CENH3_R2.fq.gz \
	    | samtools sort -@ 8 -o 04_cenh3_signal/${SP}.CENH3.raw.bam

bwa-mem2 mem -t $THREADS $GENOME \
	  02_chipseq/${SP}_Input_R1.fq.gz 02_chipseq/${SP}_Input_R2.fq.gz \
	    | samtools sort -@ 8 -o 04_cenh3_signal/${SP}.Input.raw.bam

samtools index 04_cenh3_signal/${SP}.CENH3.raw.bam
samtools index 04_cenh3_signal/${SP}.Input.raw.bam

samtools view -@ 8 -b -q 20 -F 1804 \
  04_cenh3_signal/${SP}.CENH3.raw.bam \
  > 04_cenh3_signal/${SP}.CENH3.q20.bam

samtools view -@ 8 -b -q 20 -F 1804 \
  04_cenh3_signal/${SP}.Input.raw.bam \
  > 04_cenh3_signal/${SP}.Input.q20.bam

samtools index 04_cenh3_signal/${SP}.CENH3.q20.bam
samtools index 04_cenh3_signal/${SP}.Input.q20.bam


samtools view -@ 8 -b -F 1804 \
  04_cenh3_signal/${SP}.CENH3.raw.bam \
  > 04_cenh3_signal/${SP}.CENH3.relaxed.bam

samtools view -@ 8 -b -F 1804 \
  04_cenh3_signal/${SP}.Input.raw.bam \
  > 04_cenh3_signal/${SP}.Input.relaxed.bam

samtools index 04_cenh3_signal/${SP}.CENH3.relaxed.bam
samtools index 04_cenh3_signal/${SP}.Input.relaxed.bam
