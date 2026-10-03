#!/usr/bin/env bash
set -euo pipefail

# Full sequence-based synteny workflow prototype for A17 vs R108.
# Run on lz10 after installing/activating the cent_synteny environment.
#
# Important:
# - The earlier genome_A17_vs_genome_R108.asm5.paf does not include cs tags.
# - paftools.js needs cs tags to convert PAF to MAF.
# - Therefore this script generates a second PAF with --cs=long.

source ~/.bashrc 2>/dev/null || true
conda activate cent_synteny

BASE=path/to/project/10.centromere_analysis
OUT=${BASE}/output_synthen/02_cent_synteny/genome_A17_genome_R108

cd "${OUT}"

TARGET=genome_R108
QUERY=genome_A17
THREADS=${THREADS:-24}

ln -sf ../../../output_all/00_genome/${TARGET}.fa ${TARGET}.fa
ln -sf ../../../output_all/00_genome/${QUERY}.fa ${QUERY}.fa

# UCSC mafToAxt works best when MAF source names use db.sequence naming.
# Prefix FASTA headers so MAF rows become e.g. genome_R108.Chr1 and genome_A17.Chr1.
if [ ! -s ${TARGET}.prefixed.fa ]; then
  awk -v db="${TARGET}" '/^>/{sub(/^>/, ">" db "."); print; next} {print}' ${TARGET}.fa > ${TARGET}.prefixed.fa
fi

if [ ! -s ${QUERY}.prefixed.fa ]; then
  awk -v db="${QUERY}" '/^>/{sub(/^>/, ">" db "."); print; next} {print}' ${QUERY}.fa > ${QUERY}.prefixed.fa
fi

faToTwoBit ${TARGET}.prefixed.fa ${TARGET}.2bit
faToTwoBit ${QUERY}.prefixed.fa ${QUERY}.2bit
twoBitInfo ${TARGET}.2bit ${TARGET}.sizes
twoBitInfo ${QUERY}.2bit ${QUERY}.sizes

if [ ! -s ${QUERY}_vs_${TARGET}.asm5.cs.paf ]; then
  minimap2 -x asm5 -t "${THREADS}" -c --cs=long --eqx ${TARGET}.prefixed.fa ${QUERY}.prefixed.fa > ${QUERY}_vs_${TARGET}.asm5.cs.paf
fi

paftools.js view -f maf ${QUERY}_vs_${TARGET}.asm5.cs.paf > ${QUERY}_vs_${TARGET}.asm5.raw.maf
python normalize_paftools_maf.py ${QUERY}_vs_${TARGET}.asm5.raw.maf ${QUERY}_vs_${TARGET}.asm5.maf

# UCSC mafToAxt expects a target name and query db name. Because the FASTA
# headers are prefixed as db.sequence, this should convert all chromosomes.
mafToAxt ${QUERY}_vs_${TARGET}.asm5.maf first ${QUERY} ${QUERY}_vs_${TARGET}.asm5.axt

axtSort ${QUERY}_vs_${TARGET}.asm5.axt ${QUERY}_vs_${TARGET}.asm5.sorted.axt
axtChain -linearGap=medium ${QUERY}_vs_${TARGET}.asm5.sorted.axt ${TARGET}.2bit ${QUERY}.2bit ${QUERY}_vs_${TARGET}.asm5.chain
chainSort ${QUERY}_vs_${TARGET}.asm5.chain ${QUERY}_vs_${TARGET}.asm5.sorted.chain
chainPreNet ${QUERY}_vs_${TARGET}.asm5.sorted.chain ${TARGET}.sizes ${QUERY}.sizes ${QUERY}_vs_${TARGET}.asm5.prenet.chain
chainNet ${QUERY}_vs_${TARGET}.asm5.prenet.chain ${TARGET}.sizes ${QUERY}.sizes ${TARGET}.${QUERY}.net ${QUERY}.${TARGET}.net
netSyntenic ${TARGET}.${QUERY}.net ${TARGET}.${QUERY}.syntenic.net
netSyntenic ${QUERY}.${TARGET}.net ${QUERY}.${TARGET}.syntenic.net
netToAxt ${TARGET}.${QUERY}.syntenic.net ${QUERY}_vs_${TARGET}.asm5.prenet.chain ${TARGET}.2bit ${QUERY}.2bit ${QUERY}_vs_${TARGET}.asm5.netSyntenic.axt

ls -lh \
  ${QUERY}_vs_${TARGET}.asm5.cs.paf \
  ${QUERY}_vs_${TARGET}.asm5.maf \
  ${QUERY}_vs_${TARGET}.asm5.axt \
  ${QUERY}_vs_${TARGET}.asm5.sorted.chain \
  ${TARGET}.${QUERY}.syntenic.net \
  ${QUERY}_vs_${TARGET}.asm5.netSyntenic.axt
