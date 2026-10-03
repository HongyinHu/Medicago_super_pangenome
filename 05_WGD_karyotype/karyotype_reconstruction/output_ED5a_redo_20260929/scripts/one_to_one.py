#!/usr/bin/env python3
"""Genome-wide one-to-one filter on WGDI block tables.

WGDI 0.5.1 ``-c`` only de-duplicates blocks within the same chromosome pair,
so a species region can still be claimed by blocks on several AMK chromosomes.
Here blocks are taken longest-first and kept only if at least MIN_NEW of their
genes are not yet covered, on BOTH the species side and the AMK side (the same
50 % new-coverage rule WGDI uses inside -c).
Dropped blocks are written to <out>.dropped.tsv together with the kept blocks
they conflict with, so that e.g. an AMK region present twice is visible.
Usage: one_to_one.py in.csv out.csv
"""
import sys
from collections import defaultdict

import pandas as pd

MIN_NEW = 0.5


def main(inp, out):
    df = pd.read_csv(inp)
    df["chr1"] = df["chr1"].astype(str)
    df["chr2"] = df["chr2"].astype(str)
    df = df.sort_values(["length", "pvalue"], ascending=[False, True])
    # owner[(side, chr, order)] = id of the kept block covering that gene
    owner = {}
    keep, dropped = [], []
    for idx, r in df.iterrows():
        b1 = [int(x) for x in str(r["block1"]).split("_") if x]
        b2 = [int(x) for x in str(r["block2"]).split("_") if x]
        hit1 = [owner[(1, r["chr1"], x)] for x in b1 if (1, r["chr1"], x) in owner]
        hit2 = [owner[(2, r["chr2"], x)] for x in b2 if (2, r["chr2"], x) in owner]
        new1 = 1 - len(hit1) / float(len(b1))
        new2 = 1 - len(hit2) / float(len(b2))
        if new1 >= MIN_NEW and new2 >= MIN_NEW:
            keep.append(idx)
            for x in b1:
                owner.setdefault((1, r["chr1"], x), r["id"])
            for x in b2:
                owner.setdefault((2, r["chr2"], x), r["id"])
        else:
            dropped.append({
                "id": r["id"], "chr1": r["chr1"], "start1": min(b1), "end1": max(b1),
                "chr2_AMK": r["chr2"], "start2": min(b2), "end2": max(b2),
                "length": r["length"], "ks_median": r["ks_median"], "pvalue": r["pvalue"],
                "new_cov_species": round(new1, 3), "new_cov_AMK": round(new2, 3),
                "conflict_side": ("species" if new1 < MIN_NEW else "") +
                                 ("+" if new1 < MIN_NEW and new2 < MIN_NEW else "") +
                                 ("AMK" if new2 < MIN_NEW else ""),
                "conflicting_kept_ids": ",".join(str(i) for i in sorted(set(hit1 + hit2))),
            })
    res = df.loc[keep].sort_values("id")
    res.to_csv(out, index=False)
    pd.DataFrame(dropped, columns=[
        "id", "chr1", "start1", "end1", "chr2_AMK", "start2", "end2", "length", "ks_median",
        "pvalue", "new_cov_species", "new_cov_AMK", "conflict_side", "conflicting_kept_ids",
    ]).to_csv(out + ".dropped.tsv", sep="\t", index=False)
    print("one_to_one: {} -> {} blocks".format(len(df), len(res)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
