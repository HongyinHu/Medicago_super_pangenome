#!/usr/bin/env python3
import argparse

import pandas as pd


def read_bed(path):
    return pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end"], usecols=[0, 1, 2])


def add_overlap_flag(blocks, bed, prefix):
    flag = pd.Series(False, index=blocks.index)
    chr_col = f"{prefix}_chr"
    start_col = f"{prefix}_start"
    end_col = f"{prefix}_end"
    for chrom, regions in bed.groupby("chr"):
        idx = blocks[chr_col] == chrom
        if not idx.any():
            continue
        local = pd.Series(False, index=blocks.index[idx])
        for _, region in regions.iterrows():
            local = local | (
                (blocks.loc[idx, start_col] < region["end"])
                & (blocks.loc[idx, end_col] > region["start"])
            )
        flag.loc[idx] = local
    return flag


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks", required=True)
    parser.add_argument("--query_cen", required=True)
    parser.add_argument("--target_cen", required=True)
    args = parser.parse_args()

    blocks = pd.read_csv(args.blocks, sep="\t")
    query_cen = read_bed(args.query_cen)
    target_cen = read_bed(args.target_cen)
    same_chr = blocks[blocks["query_chr"] == blocks["target_chr"]].copy()
    same_chr["query_cen_overlap"] = add_overlap_flag(same_chr, query_cen, "query")
    same_chr["target_cen_overlap"] = add_overlap_flag(same_chr, target_cen, "target")

    rows = []
    for min_len in [0, 100, 500, 1000, 2000, 5000, 10000, 20000]:
        sub = same_chr[
            (same_chr["query_len"] >= min_len)
            & (same_chr["target_len"] >= min_len)
        ].copy()
        cen = sub[sub["query_cen_overlap"] | sub["target_cen_overlap"]]
        rows.append(
            {
                "min_len": min_len,
                "same_chr_blocks": len(sub),
                "centromere_overlap_blocks": len(cen),
            }
        )

    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
