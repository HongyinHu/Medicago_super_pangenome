#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def read_fasta(path):
    records = {}
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = (header, "".join(chunks))
                header = line[1:]
                chunks = []
            elif line:
                chunks.append(line)
    if header is not None:
        records[header.split()[0]] = (header, "".join(chunks))
    return records


def top_ids(path, limit):
    ids = []
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            seq_id = line.split("\t")[1]
            if seq_id not in ids:
                ids.append(seq_id)
            if len(ids) == limit:
                break
    return ids


def write_fasta(path, species, ids, records):
    with open(path, "w") as handle:
        for seq_id in ids:
            seq = records[seq_id][1]
            handle.write(">{0}|{1}\n".format(species, seq_id))
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--blastp", required=True)
    parser.add_argument("--castor-pep", required=True)
    parser.add_argument("--r108-pep", required=True)
    parser.add_argument("--m22-pep", required=True)
    args = parser.parse_args()
    out = Path(args.outdir)
    combined = out / "top20_R108_M22_candidates.pep.fa"
    r108_ids = top_ids(out / "RcMYB106_vs_R108.blastp.tsv", 20)
    m22_ids = top_ids(out / "RcMYB106_vs_M22.blastp.tsv", 20)
    r108_records = read_fasta(args.r108_pep)
    m22_records = read_fasta(args.m22_pep)
    write_fasta(combined, "R108", r108_ids, r108_records)
    with open(combined, "a") as handle:
        for seq_id in m22_ids:
            seq = m22_records[seq_id][1]
            handle.write(">M22|{0}\n".format(seq_id))
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + "\n")
    output = out / "top20_R108_M22_candidates_vs_castor.blastp.tsv"
    subprocess.check_call([
        args.blastp, "-query", str(combined), "-subject", args.castor_pep,
        "-evalue", "1e-10", "-max_target_seqs", "10",
        "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore qcovs",
        "-out", str(output),
    ])
    seen = set()
    print("query\treciprocal_top\tpident\tqcovs\tevalue\tbitscore")
    with open(output) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields[0] in seen:
                continue
            seen.add(fields[0])
            print("\t".join([fields[0], fields[1], fields[2], fields[8], fields[6], fields[7]]))


if __name__ == "__main__":
    main()
