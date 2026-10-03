#!/usr/bin/env python3
import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--blocks", required=True)
    parser.add_argument("--min_len", type=int, default=10000)
    args = parser.parse_args()

    blocks = pd.read_csv(args.blocks, sep="\t")
    filtered = blocks[
        (blocks["query_len"] >= args.min_len)
        & (blocks["target_len"] >= args.min_len)
    ]

    print(f"all_blocks\t{len(blocks)}")
    print(f"min{args.min_len}_blocks\t{len(filtered)}")
    print("query_chr\ttarget_chr\tn")
    grouped = (
        filtered.groupby(["query_chr", "target_chr"])
        .size()
        .reset_index(name="n")
        .sort_values(["query_chr", "target_chr"])
    )
    for _, row in grouped.iterrows():
        print(f"{row['query_chr']}\t{row['target_chr']}\t{row['n']}")


if __name__ == "__main__":
    main()
