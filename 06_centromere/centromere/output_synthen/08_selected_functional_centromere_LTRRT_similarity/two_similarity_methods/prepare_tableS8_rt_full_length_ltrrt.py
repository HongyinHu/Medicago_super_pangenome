#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

import pandas as pd


GENOME_FASTA = {
    "genome_R108": "path/to/project/8.T2T_TEs_anno_new/output/genome_R108/EDTA/genome_R108.T2T.ctg.final.fa.mod",
    "genome_A17": "path/to/project/8.T2T_TEs_anno_new/output/genome_A17/EDTA/genome_A17.T2T.ctg.final.fa.mod",
    "genome_Mpo": "path/to/project/8.T2T_TEs_anno_new/output/genome_Mpo/EDTA/genome_Mpo.T2T.ctg.final.fa.mod",
    "genome_Msa": "path/to/project/6.genome_gap_fill_new/output/genome_Msa.T2T.ctg.final.fa",
    "genome_474": "path/to/project/6.genome_gap_fill_new/output/genome_474.near_T2T.ctg.final.fa",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def is_true(value):
    return str(value).strip().lower() == "true"


def wrap(seq, width=60):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def samtools_faidx(fasta, region):
    proc = subprocess.run(
        ["samtools", "faidx", fasta, region],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    seq = []
    for line in proc.stdout.splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--detail", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(args.detail, sep="\t", dtype=str)
    sub = meta[meta["in_functional_centromere"].map(is_true) & meta["has_rt"].map(is_true)].copy()
    sub["start"] = sub["start"].astype(int)
    sub["end"] = sub["end"].astype(int)
    sub["length_bp"] = sub["length_bp"].astype(int)
    sub["full_start_1based"] = sub["start"] + 1
    sub["full_end_1based"] = sub["end"]
    sub["insertion_time_mya"] = pd.to_numeric(sub["insertion_time_years"], errors="coerce") / 1_000_000.0
    sub["selected_type"] = sub["superfamily"].astype(str) + "_" + sub["clade"].astype(str)
    sub["selected_type"] = sub["selected_type"].str.replace(" ", "_", regex=False)

    records = []
    fasta_path = outdir / "TableS8_functional_centromere_RT_domain_full_length_LTRRT.fa"
    with fasta_path.open("w") as out:
        for _, row in sub.iterrows():
            genome = row["genome"]
            if genome not in GENOME_FASTA:
                raise SystemExit("No FASTA configured for {}".format(genome))
            region = "{}:{}-{}".format(row["chr"], row["full_start_1based"], row["full_end_1based"])
            seq = samtools_faidx(GENOME_FASTA[genome], region).upper()
            if row.get("strand", "") == "-":
                seq = seq.translate(COMP)[::-1]
            out.write(">{}\n".format(row["seq_id"]))
            for part in wrap(seq):
                out.write(part + "\n")
            rec = row.to_dict()
            rec["extracted_region"] = region
            rec["extracted_length_bp"] = len(seq)
            records.append(rec)

    out_meta = pd.DataFrame(records)
    out_meta.to_csv(outdir / "TableS8_functional_centromere_RT_domain_full_length_LTRRT.metadata.tsv", sep="\t", index=False)
    print("selected_records={}".format(len(out_meta)))
    print(out_meta.groupby(["label", "selected_type"]).size().reset_index(name="count").to_string(index=False))
    print("fasta={}".format(fasta_path))


if __name__ == "__main__":
    main()
