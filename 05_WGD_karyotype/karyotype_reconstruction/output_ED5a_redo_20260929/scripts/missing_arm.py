#!/usr/bin/env python3
"""Are AMK4 1-2478 (Marc) and AMK2 1-2735 (Mlup) absent from the assembled
chromosomes, or only filtered out?  Checks raw WGDI blocks, the diamond hits,
and how many annotated genes sit outside the chromosome set."""
import os

import pandas as pd

R = "path/to/project/37.karyotype_reconstruction/output_Kary"
CASES = [("Marc", "genome_Mar", "4", 1, 2478), ("Mlup", "genome_395", "2", 1, 2735)]
for lab, sub, amk, lo, hi in CASES:
    d = os.path.join(R, sub)
    conf = open(os.path.join(d, "total.conf")).read()
    bi = [l.split("=", 1)[1].strip() for l in conf.splitlines() if l.strip().startswith("savefile") and "block_information" in l][0]
    raw = pd.read_csv(os.path.join(d, bi))
    raw = raw[raw.chr2.astype(str) == amk]
    cover = []
    for _, r in raw.iterrows():
        b2 = [int(x) for x in r.block2.split("_")]
        n = sum(lo <= x <= hi for x in b2)
        if n >= 5:
            cover.append((r.chr1, r.length, n, r.ks_median))
    print("##", lab, "AMK{}:{}-{} -> raw blocks with >=5 genes there:".format(amk, lo, hi))
    for c in sorted(cover, key=lambda x: -x[2])[:8]:
        print("   chr{} len{} genes_in_region{} ks{:.3f}".format(*c))
    lensf = [l.split("=", 1)[1].strip() for l in conf.splitlines() if l.strip().startswith("lens1")][0]
    print("   lens1:", open(lensf if os.path.isabs(lensf) else os.path.join(d, lensf)).read().replace("\n", " | "))
    # diamond best hits of AMK genes in the region -> where are they in the species?
    gff2 = pd.read_csv(os.path.join(d, "aak.gff"), sep="\t", header=None)
    ids = set(gff2[(gff2[0].astype(str) == amk) & (gff2[5] >= lo) & (gff2[5] <= hi)][1])
    bl = pd.read_csv(os.path.join(d, "all.blastp.txt"), sep="\t", header=None, usecols=[0, 1, 11])
    bl = bl[bl[1].isin(ids)].sort_values(11, ascending=False).drop_duplicates(1)
    gffp = [l.split("=", 1)[1].strip() for l in conf.splitlines() if l.strip().startswith("gff1")][0]
    g1 = pd.read_csv(gffp if os.path.isabs(gffp) else os.path.join(d, gffp), sep="\t", header=None, usecols=[0, 1])
    g1.columns = ["chr", "gene"]
    hit = bl.merge(g1, left_on=0, right_on="gene", how="left")
    print("   AMK genes in region:", len(ids), " with any hit:", len(bl),
          " best-hit chromosome counts:", hit.chr.value_counts().head(8).to_dict())
