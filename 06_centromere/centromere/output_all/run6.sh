SP=genome_474

# for SP in sp1 sp2 sp3 sp4 sp5
# do
cat > ${SP}.summary_to_bed.py <<'PY'
import pandas as pd
import sys

sp = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

df = pd.read_csv(infile)
df = df.loc[:, ~df.columns.str.contains("^Unnamed|^$")]

records = []

for i, row in df.iterrows():
    chrom = str(row["name"])
    start = int(row["start"])
    end = int(row["end"]) + 1

    monomer_len = int(row["most.freq.value.N"])
    array_len = end - start

    rep_class = str(row["class"]) if "class" in df.columns else "NA"
    if rep_class == "NA" or rep_class == "nan" or rep_class == "":
        repeat_family = f"{sp}_TRASH_{monomer_len}bp_{i+1:06d}"
    else:
        repeat_family = rep_class

    records.append([
        chrom,
        start,
        end,
        repeat_family,
        monomer_len,
        array_len
    ])

pd.DataFrame(records).to_csv(outfile, sep="\t", header=False, index=False)
PY

python ${SP}.summary_to_bed.py \
  ${SP} \
  03_repeat/${SP}.TRASH/Summary.of.repetitive.regions.${SP}.fa.csv \
  03_repeat/${SP}.TRASH/${SP}.TRASH_arrays.bed 

sort -k1,1 -k2,2n 03_repeat/${SP}.TRASH/${SP}.TRASH_arrays.bed > 03_repeat/${SP}.TRASH/${SP}.TRASH_arrays.sorted.bed


# 假设你已经把 TRASH summary 转成 BED
# 格式：chr start end repeat_family monomer_len array_len

bedtools intersect \
  -a 03_repeat/${SP}.TRASH/${SP}.TRASH_arrays.sorted.bed \
  -b 04_cenh3_signal/${SP}.macs2/${SP}.CENH3.q20_peaks.broadPeak \
  -wa -wb \
  > 05_centromere_define/${SP}.TRASH_overlap_CENH3peak.tsv
