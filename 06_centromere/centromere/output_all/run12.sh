#提取主要单体的bed文件
SP=genome_474
#monomer=Sat48

#for monomer in Sat344_345_346 Sat691_692 Sat138_139 Sat51 Sat170 Sat288

#for monomer in  CenSat583_584 Sat187_188 Sat143_144_145 Sat168 Sat60 Sat48 Sat185 Sat185_187_188 Sat187 Sat188

for monomer in Sat179 Sat49 Sat30 Sat54 Sat48 Sat45 Sat50 Sat29 Sat317 Sat104 Sat345
do 
  echo "now process ${monomer}======================================================================="

  GENOME=00_genome/${SP}.fa
  CENSAT_FA=05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.final_consensus.fa
  CENH3_DOMAIN=04_cenh3_signal/${SP}.macs2/${SP}.CENH3.q20_peaks.broadPeak.bed


  ls -lh ${GENOME}
  ls -lh ${CENSAT_FA}
  seqkit stat ${CENSAT_FA}

  mkdir -p 05_centromere_define/${SP}.satellite_bed

  makeblastdb \
    -in ${GENOME} \
    -dbtype nucl \
    -out 05_centromere_define/${SP}.satellite_bed/${SP}.blastdb

  blastn \
    -query ${CENSAT_FA} \
    -db 05_centromere_define/${SP}.satellite_bed/${SP}.blastdb \
    -out 05_centromere_define/${SP}.satellite_bed/${SP}_${monomer}_vs_genome.blastn.tsv \
    -outfmt "6 qseqid sseqid pident length qlen sstart send evalue bitscore" \
    -num_threads 16 \
    -task blastn-short \
    -word_size 7 \
    -dust no \
    -soft_masking false


  # 回比筛选条件
  # identity >= 95%
  # query coverage >= 95%

  awk -v monomer="$monomer" 'BEGIN{OFS="\t"}
  {
      cov=$4/$5;
      if ($3>=70 && cov>=0.50) {
          start=($6<$7?$6:$7)-1;
          end=($6>$7?$6:$7);
          print $2,start,end,monomer;
      }
  }' 05_centromere_define/${SP}.satellite_bed/${SP}_${monomer}_vs_genome.blastn.tsv \
  | sort -k1,1 -k2,2n \
  > 05_centromere_define/${SP}.satellite_bed/${SP}_${monomer}.80_80.bed

  bedtools merge \
    -i 05_centromere_define/${SP}.satellite_bed/${SP}_${monomer}.80_80.bed \
    -d 1000 \
    -c 4 -o distinct \
  > 05_centromere_define/${SP}.satellite_bed/${SP}_${monomer}.80_80.merge1k.bed

  # 只想画功能着丝粒附近的 CenSat168_169，可以和 CENH3 domain 或 MACS2 peak 取交集
  bedtools intersect \
    -a 05_centromere_define/${SP}.satellite_bed/${SP}_${monomer}.80_80.merge1k.bed \
    -b ${CENH3_DOMAIN} \
    -wa -u \
  >  05_centromere_define/${SP}.satellite_bed/${SP}_${monomer}.80_80.merge1k.CENH3_overlap.bed
done