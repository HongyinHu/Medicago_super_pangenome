#!/bin/bash
set -euo pipefail

SP=genome_474

GENOME=00_genome/${SP}.fa
TRASH_GFF=03_repeat/${SP}.TRASH/TRASH_${SP}.fa.gff
ARRAY_SUMMARY=05_centromere_define/${SP}.TRASH_CENH3_array_summary.tsv
GROUPS_file=05_centromere_define/${SP}.candidate_groups.tsv
PREFIX=${SP}

TOP_N_ARRAYS=3
SAMPLE_N=300
THREADS=16

OUTDIR=05_centromere_define/${PREFIX}.satellite_validation
mkdir -p ${OUTDIR}

echo "[0] Check input files"

for f in ${GENOME} ${TRASH_GFF} ${ARRAY_SUMMARY} ${GROUPS_file}
do
    if [ ! -s "$f" ]; then
        echo "ERROR: missing file: $f"
        exit 1
    fi
done

echo "[0.1] Prepare candidate groups"

if [ ! -s ${GROUPS_file} ]; then
    echo "ERROR: missing file: ${GROUPS_file}"
    exit 1
fi


echo "[1] Index genome"
samtools faidx ${GENOME}

echo "[2] Convert TRASH GFF to monomer BED"

awk 'BEGIN{OFS="\t"}
$0 !~ /^#/ && NF >= 9 {
    id=$1":"$4"-"$5"|"NR;
    print $1, $4-1, $5, id, ".", $7
}' ${TRASH_GFF} \
| sort -k1,1 -k2,2n \
> ${OUTDIR}/${PREFIX}.TRASH_monomers.bed

MONOMER_BED=${OUTDIR}/${PREFIX}.TRASH_monomers.bed

echo "[3] Select candidate arrays for each monomer group"

export ARRAY_SUMMARY
export GROUPS_file
export OUTDIR
export TOP_N_ARRAYS

python - <<'PY'
import os
import pandas as pd

array_summary = os.environ["ARRAY_SUMMARY"]
groups_file = os.environ["GROUPS_file"]
outdir = os.environ["OUTDIR"]
top_n = int(os.environ["TOP_N_ARRAYS"])

df = pd.read_csv(array_summary, sep="\t")
groups = pd.read_csv(groups_file, sep="\t")

required_cols = [
    "TRASH_chr", "TRASH_start", "TRASH_end", "TRASH_id",
    "monomer_len", "array_len", "CENH3_overlap_bp",
    "copy_number_est", "CENH3_overlap_ratio"
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise SystemExit(f"Missing columns in array summary: {missing}")

records = []

for _, g in groups.iterrows():
    group = str(g["group"])
    lens = [int(x) for x in str(g["monomer_lens"]).split(",")]

    sub = df[df["monomer_len"].isin(lens)].copy()

    if sub.empty:
        print(f"[WARN] no arrays found for {group}: {lens}")
        continue

    sub = sub.sort_values(
        ["CENH3_overlap_bp", "array_len", "copy_number_est"],
        ascending=False
    ).head(top_n)

    group_dir = os.path.join(outdir, group)
    os.makedirs(group_dir, exist_ok=True)

    bed_file = os.path.join(group_dir, f"{group}.candidate_arrays.bed")
    info_file = os.path.join(group_dir, f"{group}.candidate_arrays.tsv")

    sub[[
        "TRASH_chr", "TRASH_start", "TRASH_end", "TRASH_id",
        "monomer_len", "array_len", "CENH3_overlap_bp",
        "copy_number_est", "CENH3_overlap_ratio"
    ]].to_csv(info_file, sep="\t", index=False)

    with open(bed_file, "w") as out:
        for _, r in sub.iterrows():
            out.write(
                f"{r['TRASH_chr']}\t{int(r['TRASH_start'])}\t{int(r['TRASH_end'])}\t"
                f"{r['TRASH_id']}\t{int(r['monomer_len'])}\t{int(r['array_len'])}\n"
            )

    for _, r in sub.iterrows():
        records.append([
            group,
            r["TRASH_chr"],
            int(r["TRASH_start"]),
            int(r["TRASH_end"]),
            r["TRASH_id"],
            int(r["monomer_len"]),
            int(r["array_len"]),
            r["CENH3_overlap_bp"],
            r["copy_number_est"],
            r["CENH3_overlap_ratio"]
        ])

allout = pd.DataFrame(records, columns=[
    "group", "chr", "start", "end", "TRASH_id", "monomer_len",
    "array_len", "CENH3_overlap_bp", "copy_number_est",
    "CENH3_overlap_ratio"
])

allout.to_csv(
    os.path.join(outdir, "all_selected_candidate_arrays.tsv"),
    sep="\t",
    index=False
)

print(f"selected arrays written to {outdir}/all_selected_candidate_arrays.tsv")
PY

echo "[4] Build genome BLAST database"

makeblastdb \
  -in ${GENOME} \
  -dbtype nucl \
  -out ${OUTDIR}/${PREFIX}.blastdb \
  > ${OUTDIR}/${PREFIX}.makeblastdb.log 2>&1

echo "[5] Prepare consensus script"

cat > ${OUTDIR}/make_consensus_from_alignment.py <<'PY'
import sys
from collections import Counter

aln = sys.argv[1]
out = sys.argv[2]
name = sys.argv[3]

seqs = []
seq = []

with open(aln) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq:
                seqs.append("".join(seq).upper())
            seq = []
        else:
            seq.append(line)
    if seq:
        seqs.append("".join(seq).upper())

if not seqs:
    raise SystemExit("No sequences found in alignment")

length = max(len(s) for s in seqs)
consensus = []

for i in range(length):
    bases = []
    for s in seqs:
        if i < len(s):
            b = s[i]
            if b in "ACGT":
                bases.append(b)

    if not bases:
        continue

    c = Counter(bases)
    base, count = c.most_common(1)[0]

    if count / len(bases) >= 0.5:
        consensus.append(base)
    else:
        consensus.append("N")

consensus = "".join(consensus).replace("-", "")

with open(out, "w") as f:
    f.write(f">{name}|len={len(consensus)}|ncopy={len(seqs)}\n")
    for i in range(0, len(consensus), 80):
        f.write(consensus[i:i+80] + "\n")

print(f"{name}\tconsensus_len={len(consensus)}\tncopy={len(seqs)}")
PY

echo "[6] Process selected candidate arrays"

tail -n +2 ${OUTDIR}/all_selected_candidate_arrays.tsv | while IFS=$'\t' read group chr start end trash_id mlen array_len cenh3_overlap copy_est overlap_ratio
do
    safe_id=$(echo ${trash_id} | sed 's/[^A-Za-z0-9_.-]/_/g')
    array_name=${group}.${safe_id}.${chr}_${start}_${end}
    work=${OUTDIR}/${group}/${array_name}
    mkdir -p ${work}

    echo "Processing ${array_name}"

    echo -e "${chr}\t${start}\t${end}\t${array_name}" \
      > ${work}/${array_name}.array.bed

    bedtools intersect \
      -a ${MONOMER_BED} \
      -b ${work}/${array_name}.array.bed \
      -wa \
      > ${work}/${array_name}.monomers.raw.bed

    awk -v m=${mlen} 'BEGIN{OFS="\t"}
    {
        len=$3-$2;
        if (len >= m*0.70 && len <= m*1.30) print $0
    }' ${work}/${array_name}.monomers.raw.bed \
    > ${work}/${array_name}.monomers.filtered.bed

    n_monomer=$(wc -l < ${work}/${array_name}.monomers.filtered.bed)

    if [ "${n_monomer}" -lt 5 ]; then
        echo "[WARN] too few monomer copies for ${array_name}: ${n_monomer}"
        continue
    fi

    bedtools getfasta \
      -fi ${GENOME} \
      -bed ${work}/${array_name}.monomers.filtered.bed \
      -s \
      -name \
      -fo ${work}/${array_name}.monomers.fa

    n_fa=$(grep -c "^>" ${work}/${array_name}.monomers.fa || true)

    if [ "${n_fa}" -gt "${SAMPLE_N}" ]; then
        seqkit sample \
          -n ${SAMPLE_N} \
          ${work}/${array_name}.monomers.fa \
          > ${work}/${array_name}.monomers.sample.fa
    else
        cp ${work}/${array_name}.monomers.fa \
          ${work}/${array_name}.monomers.sample.fa
    fi

    mafft \
      --auto \
      --adjustdirectionaccurately \
      --thread ${THREADS} \
      ${work}/${array_name}.monomers.sample.fa \
      > ${work}/${array_name}.monomers.sample.aln.fa \
      2> ${work}/${array_name}.mafft.log

    python ${OUTDIR}/make_consensus_from_alignment.py \
      ${work}/${array_name}.monomers.sample.aln.fa \
      ${work}/${array_name}.consensus.fa \
      ${array_name} \
      > ${work}/${array_name}.consensus.stat.txt

    makeblastdb \
      -in ${work}/${array_name}.consensus.fa \
      -dbtype nucl \
      -out ${work}/${array_name}.consensus.db \
      > ${work}/${array_name}.consensus.makeblastdb.log 2>&1

    blastn \
      -query ${work}/${array_name}.consensus.fa \
      -db ${work}/${array_name}.consensus.db \
      -task blastn-short \
      -dust no \
      -word_size 7 \
      -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
      -out ${work}/${array_name}.selfblast.tsv \
      -num_threads ${THREADS}

    if command -v trf >/dev/null 2>&1; then
        cd ${work}
        trf ${array_name}.consensus.fa 2 7 7 80 10 50 2000 -d -h \
          > ${array_name}.trf.log 2>&1 || true
        cd - >/dev/null
    else
        echo "[WARN] TRF not found, skip TRF for ${array_name}"
    fi

    blastn \
      -query ${work}/${array_name}.consensus.fa \
      -db ${OUTDIR}/${PREFIX}.blastdb \
      -out ${work}/${array_name}.blast_genome.tsv \
      -outfmt "6 qseqid sseqid pident length qlen sstart send evalue bitscore" \
      -num_threads ${THREADS}

    awk 'BEGIN{OFS="\t"}
    {
        cov=$4/$5;
        if ($3>=95 && cov>=0.95) {
            start=($6<$7?$6:$7)-1;
            end=($6>$7?$6:$7);
            print $2,start,end,$1,$3,cov;
        }
    }' ${work}/${array_name}.blast_genome.tsv \
    | sort -k1,1 -k2,2n \
    > ${work}/${array_name}.blast_genome.95_95.bed

    awk 'BEGIN{OFS="\t"}
    {
        cov=$4/$5;
        if ($3>=80 && cov>=0.80) {
            start=($6<$7?$6:$7)-1;
            end=($6>$7?$6:$7);
            print $2,start,end,$1,$3,cov;
        }
    }' ${work}/${array_name}.blast_genome.tsv \
    | sort -k1,1 -k2,2n \
    > ${work}/${array_name}.blast_genome.80_80.bed

    c95=$(wc -l < ${work}/${array_name}.blast_genome.95_95.bed)
    c80=$(wc -l < ${work}/${array_name}.blast_genome.80_80.bed)

    echo -e "array_name\tgroup\tchr\tstart\tend\tmonomer_len\tarray_len\textracted_copy\tconsensus_fa\tblast95_count\tblast80_count" \
      > ${work}/${array_name}.final_stat.tsv

    echo -e "${array_name}\t${group}\t${chr}\t${start}\t${end}\t${mlen}\t${array_len}\t${n_fa}\t${work}/${array_name}.consensus.fa\t${c95}\t${c80}" \
      >> ${work}/${array_name}.final_stat.tsv

done

echo "[7] Collect final stats"

echo -e "array_name\tgroup\tchr\tstart\tend\tmonomer_len\tarray_len\textracted_copy\tconsensus_fa\tblast95_count\tblast80_count" \
> ${OUTDIR}/all_candidate_consensus.final_stat.tsv

find ${OUTDIR} -name "*.final_stat.tsv" | sort | while read f
do
    tail -n +2 "$f"
done >> ${OUTDIR}/all_candidate_consensus.final_stat.tsv

echo "[8] Collect all consensus FASTA"

find ${OUTDIR} -name "*.consensus.fa" | sort | while read f
do
    cat "$f"
done > ${OUTDIR}/${PREFIX}.all_candidate_consensus.fa

echo "Done."
echo "Main outputs:"
echo "${OUTDIR}/all_selected_candidate_arrays.tsv"
echo "${OUTDIR}/all_candidate_consensus.final_stat.tsv"
echo "${OUTDIR}/${PREFIX}.all_candidate_consensus.fa"
