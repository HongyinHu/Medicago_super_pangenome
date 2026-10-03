#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd


SELECTED_BY_LABEL = {
    "R108": {("Gypsy", "CRM"), ("Gypsy", "Tekay")},
    "A17": {("Gypsy", "CRM"), ("Gypsy", "Tekay")},
    "Mpo": {("Copia", "SIRE"), ("Gypsy", "Athila")},
    "Msa": {("Copia", "SIRE"), ("Gypsy", "CRM")},
    "474": {("Copia", "SIRE"), ("Gypsy", "Athila"), ("Gypsy", "Ogre"), ("Gypsy", "Tekay")},
}


def truthy(value):
    return str(value).strip().lower() in {"true", "1", "yes"}


def rename_fasta(raw_fasta, renamed_fasta, seq_ids):
    idx = -1
    with open(raw_fasta) as inp, open(renamed_fasta, "w") as out:
        for line in inp:
            if line.startswith(">"):
                idx += 1
                if idx >= len(seq_ids):
                    raise RuntimeError(f"More FASTA records than seq IDs in {raw_fasta}")
                out.write(f">{seq_ids[idx]}\n")
            else:
                out.write(line)
    if idx + 1 != len(seq_ids):
        raise RuntimeError(f"FASTA records ({idx + 1}) != seq IDs ({len(seq_ids)}) in {raw_fasta}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--detail")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--rename-raw-fasta")
    parser.add_argument("--renamed-fasta")
    parser.add_argument("--seq-ids")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.rename_raw_fasta:
        if not args.renamed_fasta or not args.seq_ids:
            raise SystemExit("--renamed-fasta and --seq-ids are required with --rename-raw-fasta")
        seq_ids = [line.strip().split()[0] for line in open(args.seq_ids) if line.strip()]
        rename_fasta(args.rename_raw_fasta, args.renamed_fasta, seq_ids)
        return

    if not args.detail:
        raise SystemExit("--detail is required unless --rename-raw-fasta is used")

    df = pd.read_csv(args.detail, sep="\t")
    df["label"] = df["label"].astype(str).str.strip()
    df["superfamily"] = df["superfamily"].astype(str).str.strip()
    df["clade"] = df["clade"].astype(str).str.strip()
    df = df[df["in_functional_centromere"].map(truthy)].copy()
    df = df[
        df.apply(
            lambda r: (str(r["superfamily"]), str(r["clade"])) in SELECTED_BY_LABEL.get(str(r["label"]), set()),
            axis=1,
        )
    ].copy()

    df["selected_type"] = df["superfamily"].astype(str) + "_" + df["clade"].astype(str)
    df["full_start_1based"] = df["start"].astype(int) + 1
    df["full_end_1based"] = df["end"].astype(int)
    df["faidx_region"] = df["chr"].astype(str) + ":" + df["full_start_1based"].astype(str) + "-" + df["full_end_1based"].astype(str)
    df["insertion_time_mya"] = pd.to_numeric(df["insertion_time_years"], errors="coerce") / 1_000_000

    counts = {}
    seq_ids = []
    for _, row in df.iterrows():
        label = str(row["label"])
        typ = str(row["selected_type"])
        counts[(label, typ)] = counts.get((label, typ), 0) + 1
        safe_type = re.sub(r"[^A-Za-z0-9_.-]+", "_", typ)
        seq_ids.append(
            f"{label}__{safe_type}__{row['chr']}__{int(row['full_start_1based'])}_{int(row['full_end_1based'])}__N{counts[(label, typ)]}"
        )
    df["seq_id"] = seq_ids

    meta_path = outdir / "selected_functional_centromere_full_length_LTRRT.metadata.tsv"
    df.to_csv(meta_path, sep="\t", index=False)

    for label, sub in df.groupby("label", sort=False):
        subdir = outdir / str(label)
        subdir.mkdir(parents=True, exist_ok=True)
        sub[["faidx_region"]].to_csv(subdir / f"{label}.full_length_LTRRT.regions.txt", sep="\t", index=False, header=False)
        sub[["seq_id"]].to_csv(subdir / f"{label}.full_length_LTRRT.seq_ids.txt", sep="\t", index=False, header=False)

    summary = (
        df.groupby(["label", "superfamily", "clade", "selected_type"])
        .size()
        .reset_index(name="count")
        .sort_values(["selected_type", "label"])
    )
    summary.to_csv(outdir / "selected_functional_centromere_full_length_LTRRT.counts.tsv", sep="\t", index=False)

    print(f"metadata\t{meta_path}")
    print(f"total_selected\t{len(df)}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
