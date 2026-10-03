#添加trf repeat 轨道

#!/bin/bash
set -euo pipefail

SP=genome_474

GENOME=00_genome/${SP}.fa
OUTDIR=05_centromere_define/${SP}.TRF

THREADS=30

mkdir -p ${OUTDIR}/chr_fa
mkdir -p ${OUTDIR}/trf_raw
mkdir -p ${OUTDIR}/bed

echo "[1] Index genome"
samtools faidx ${GENOME}

cut -f1,2 ${GENOME}.fai > ${OUTDIR}/${SP}.chrom.sizes

echo "[2] Extract each chromosome"

while read chr len
do
    echo "Extract ${chr}"
    samtools faidx ${GENOME} ${chr} > ${OUTDIR}/chr_fa/${chr}.fa
done < ${OUTDIR}/${SP}.chrom.sizes

echo "[3] Run TRF for each chromosome"

while read chr len
do
    echo "Run TRF on ${chr}"
    cd ${OUTDIR}/trf_raw

    trf ../chr_fa/${chr}.fa 2 7 7 80 10 50 2000 -d -h

    cd -
done < ${OUTDIR}/${SP}.chrom.sizes

echo "[4] Parse TRF dat to BED"

cat > ${OUTDIR}/parse_trf_dat_to_bed.py << 'PY'
import sys
import glob
import os

out = sys.argv[1]

with open(out, "w") as w:
    for dat in glob.glob("trf_raw/*.dat"):
        chrom = None

        with open(dat) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith("Sequence:"):
                    chrom = line.split(":", 1)[1].strip().split()[0]
                    continue

                parts = line.split()
                if chrom is None:
                    continue

                if len(parts) < 14:
                    continue

                if not parts[0].isdigit():
                    continue

                start = int(parts[0])
                end = int(parts[1])
                period = parts[2]
                copy_number = parts[3]
                consensus_size = parts[4]
                perc_match = parts[5]
                perc_indel = parts[6]
                score = parts[7]
                entropy = parts[12]

                bed_start = start - 1
                bed_end = end
                array_len = end - start + 1

                name = f"TRF_{period}bp_{chrom}_{start}_{end}"

                w.write(
                    f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t{score}\t.\t"
                    f"{period}\t{copy_number}\t{consensus_size}\t"
                    f"{perc_match}\t{perc_indel}\t{array_len}\t{entropy}\n"
                )
PY

cd ${OUTDIR}

python parse_trf_dat_to_bed.py bed/${SP}.TRF.all.bed

sort -k1,1 -k2,2n bed/${SP}.TRF.all.bed > bed/${SP}.TRF.all.sorted.bed

echo "[5] Make filtered BED files for plotting"

# 所有长度 >=100 bp 且 TRF score >=50 的重复
awk 'BEGIN{OFS="\t"} $12>=100 && $5>=50 {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13
}' bed/${SP}.TRF.all.sorted.bed \
> bed/${SP}.TRF.len100.score50.bed

# 整条染色体图推荐更严格一点，只画长度 >=500 bp 的 TRF repeat
awk 'BEGIN{OFS="\t"} $12>=500 && $5>=100 {
    print $1,$2,$3,$4
}' bed/${SP}.TRF.all.sorted.bed \
> bed/${SP}.TRF.len500.score100.plot.bed

# 局部放大图可以用长度 >=100 bp 的版本
awk 'BEGIN{OFS="\t"} $12>=100 && $5>=50 {
    print $1,$2,$3,$4
}' bed/${SP}.TRF.all.sorted.bed \
> bed/${SP}.TRF.len100.score50.plot.bed

echo "[6] Make TRF density bedGraph and bigWig"

bedtools makewindows \
  -g ${SP}.chrom.sizes \
  -w 100000 \
> bed/${SP}.100kb.windows.bed

bedtools coverage \
  -a bed/${SP}.100kb.windows.bed \
  -b bed/${SP}.TRF.len100.score50.bed \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' \
> bed/${SP}.TRF.coverage.100kb.bedGraph

bedGraphToBigWig \
  bed/${SP}.TRF.coverage.100kb.bedGraph \
  ${SP}.chrom.sizes \
  bed/${SP}.TRF.coverage.100kb.bw

echo "Done."
echo "Main output files:"
echo "${OUTDIR}/bed/${SP}.TRF.len500.score100.plot.bed"
echo "${OUTDIR}/bed/${SP}.TRF.len100.score50.plot.bed"
echo "${OUTDIR}/bed/${SP}.TRF.coverage.100kb.bw"