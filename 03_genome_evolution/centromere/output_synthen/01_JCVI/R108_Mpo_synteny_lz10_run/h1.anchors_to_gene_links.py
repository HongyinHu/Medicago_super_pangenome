#!/usr/bin/env python3
import argparse
import pandas as pd


def read_coord(path, prefix):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            f"{prefix}_gene",
            f"{prefix}_chr",
            f"{prefix}_start",
            f"{prefix}_end",
            f"{prefix}_strand",
        ],
    )


def read_anchors(path):
    rows = []
    block_id = 0

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                block_id += 1
                continue

            parts = line.split()
            if len(parts) >= 2:
                rows.append([parts[0], parts[1], f"block{block_id}"])

    return pd.DataFrame(rows, columns=["sa_gene", "sm_gene", "block_id"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchors", required=True)
    parser.add_argument("--sa_coord", required=True)
    parser.add_argument("--sm_coord", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    sa = read_coord(args.sa_coord, "sa")
    sm = read_coord(args.sm_coord, "sm")
    pairs = read_anchors(args.anchors)

    df = pairs.merge(sa, on="sa_gene", how="inner")
    df = df.merge(sm, on="sm_gene", how="inner")

    out = df[
        [
            "sa_chr",
            "sa_start",
            "sa_end",
            "sm_chr",
            "sm_start",
            "sm_end",
            "block_id",
            "sa_gene",
            "sm_gene",
        ]
    ]

    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()