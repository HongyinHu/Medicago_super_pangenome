#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

import pandas as pd


def parse_locus(value):
    m = re.match(r"^([^:]+):(\d+)\.\.(\d+)$", value)
    if not m:
        return None
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    return chrom, min(start, end), max(start, end)


def parse_internal(value):
    m = re.match(r"^IN:(\d+)\.\.(\d+)$", value)
    if not m:
        return None
    start = int(m.group(1))
    end = int(m.group(2))
    return min(start, end), max(start, end)


def read_passlist(path, label):
    rows = []
    row_num = 0
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 10:
                continue
            locus = parse_locus(f[0])
            internal = parse_internal(f[6]) if len(f) > 6 else None
            if locus is None or internal is None:
                continue
            chrom, ltr_start, ltr_end = locus
            if not re.match(r"^Chr[0-9]+$", chrom):
                continue
            internal_start, internal_end = internal
            row_num += 1
            seq_id = f"{label}__{chrom}__LTR{ltr_start}_{ltr_end}__INT{internal_start}_{internal_end}__N{row_num}"
            rows.append({
                "seq_id": seq_id,
                "faidx_region": f"{chrom}:{internal_start}-{internal_end}",
                "chr": chrom,
                "start": ltr_start - 1,
                "end": ltr_end,
                "length_bp": ltr_end - ltr_start + 1,
                "ltr_locus": f[0],
                "internal_start": internal_start,
                "internal_end": internal_end,
                "internal_length_bp": internal_end - internal_start + 1,
                "pass_category": f[1] if len(f) > 1 else "NA",
                "motif": f[2] if len(f) > 2 else "NA",
                "tsd": f[3] if len(f) > 3 else "NA",
                "strand": f[8] if len(f) > 8 else "NA",
                "pass_superfamily": f[9] if len(f) > 9 else "unknown",
                "pass_TE_type": f[10] if len(f) > 10 else "NA",
                "ltr_identity": f[7] if len(f) > 7 else "NA",
                "insertion_time_years": f[11] if len(f) > 11 else "NA",
            })
    return rows


def rename_fasta(raw_fasta, renamed_fasta, seq_ids):
    idx = -1
    with open(raw_fasta) as inp, open(renamed_fasta, "w") as out:
        for line in inp:
            if line.startswith(">"):
                idx += 1
                if idx >= len(seq_ids):
                    raise RuntimeError(f"More FASTA records than metadata rows in {raw_fasta}")
                out.write(f">{seq_ids[idx]}\n")
            else:
                out.write(line)
    if idx + 1 != len(seq_ids):
        raise RuntimeError(f"FASTA records ({idx + 1}) != metadata rows ({len(seq_ids)}) in {raw_fasta}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--rename-raw-fasta")
    parser.add_argument("--renamed-fasta")
    parser.add_argument("--meta")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.rename_raw_fasta:
        if not args.meta or not args.renamed_fasta:
            raise SystemExit("--meta and --renamed-fasta are required with --rename-raw-fasta")
        meta = pd.read_csv(args.meta, sep="\t")
        rename_fasta(args.rename_raw_fasta, args.renamed_fasta, meta["seq_id"].tolist())
        return

    if not args.config:
        raise SystemExit("--config is required unless --rename-raw-fasta is used")

    config = pd.read_csv(args.config, sep="\t")
    for _, cfg in config.iterrows():
        label = cfg["label"]
        subdir = outdir / label
        subdir.mkdir(parents=True, exist_ok=True)
        rows = read_passlist(cfg["pass_list"], label)
        meta = subdir / f"{label}.passlist_internal.meta.tsv"
        regions = subdir / f"{label}.passlist_internal.regions.txt"
        pd.DataFrame(rows).to_csv(meta, sep="\t", index=False)
        with open(regions, "w") as handle:
            for row in rows:
                handle.write(row["faidx_region"] + "\n")
        print(f"{label}\t{len(rows)}\t{meta}\t{regions}")


if __name__ == "__main__":
    main()
