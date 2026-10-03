#指定候选单体生成一致性序列
#执行环境CENH3_env

SP=genome_474
#monomer=Sat50

#for monomer in Sat344_345_346 Sat691_692 Sat138_139 Sat51 Sat170 Sat288
#for monomer in CenSat583_584 Sat187_188 Sat143_144_145 Sat168 Sat60 Sat48 Sat185 Sat185_187_188 Sat187 Sat188
for monomer in Sat179 Sat49 Sat30 Sat54 Sat48 Sat45 Sat50 Sat29 Sat317 Sat104 Sat345
do 

  find 05_centromere_define/${SP}.satellite_validation/${monomer} \
    -name "*.monomers.fa" \
    -exec cat {} \; \
    > 05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.all_monomers.fa

  seqkit sample -n 1000 \
    05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.all_monomers.fa \
    > 05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.sample1000.fa


  mafft --auto --adjustdirectionaccurately --thread 16 \
    05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.sample1000.fa \
    > 05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.sample1000.aln.fa

  python 05_centromere_define/${SP}.satellite_validation/make_consensus_from_alignment.py \
    05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.sample1000.aln.fa \
    05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.final_consensus.fa \
    ${SP}_${monomer}

  # trf 05_centromere_define/${SP}.satellite_validation/${monomer}/${SP}_${monomer}.final_consensus.fa \
  #   2 7 7 80 10 50 2000 -d -h
done