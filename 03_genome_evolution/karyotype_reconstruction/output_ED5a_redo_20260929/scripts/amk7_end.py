#!/usr/bin/env python3
"""Where does the distal end of AMK7 (genes 6200-6466) sit in each genome?
Uses all raw WGDI blocks with ortholog-level Ks, >= 10 genes."""
import os

import pandas as pd

OUT = ("path/to/project/37.karyotype_reconstruction/"
       "output_ED5a_redo_20260929")
W = os.path.join(OUT, "04_wgdi_p0.2")
man = pd.read_csv(os.path.join(W, "run_manifest.tsv"), sep="\t").set_index("label")
sp = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t").sort_values("order")
for _, s in sp.iterrows():
    raw = pd.read_csv(os.path.join(W, s.label, "in.blockinfo_raw"))
    raw = raw[(raw.chr2.astype(str) == "7") & (raw.length >= 10) & (raw.ks_median <= man.loc[s.label, "ks_cut"])]
    hits = []
    for _, r in raw.iterrows():
        b2 = list(map(int, r.block2.split("_")))
        b1 = list(map(int, r.block1.split("_")))
        lo, hi = min(b2), max(b2)
        if hi >= 6200:
            hits.append("chr{}:{}-{} <- AMK7:{}-{} ({}g, Ks {:.3f})".format(
                r.chr1, min(b1), max(b1), lo, hi, r.length, r.ks_median))
    print("{:9s} {}".format(s.label, " | ".join(hits) if hits else "-"))
