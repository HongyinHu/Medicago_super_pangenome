#!/usr/bin/env python3
"""Does the chr7 block of Mrut/Marc still contain AMK7 genes 6343-6406
(i.e. copy kept on chr7 = duplication) or skip them (moved = translocation)?"""
import os

import pandas as pd

W = ("path/to/project/37.karyotype_reconstruction/"
     "output_ED5a_redo_20260929/04_wgdi_p0.2")
for lab in ("Mrut", "Marc", "Mlan", "Medg", "Mtru_R108"):
    raw = pd.read_csv(os.path.join(W, lab, "in.blockinfo_raw"))
    raw = raw[(raw.chr2.astype(str) == "7") & (raw.length >= 10) & (raw.ks_median < 0.3)]
    for _, r in raw.iterrows():
        b2 = list(map(int, r.block2.split("_")))
        inside = [x for x in b2 if 6343 <= x <= 6406]
        near = [x for x in b2 if 6300 <= x <= 6466]
        if near:
            print("{:9s} chr{} block{} ({}g): AMK7 6343-6406 genes in block = {}, 6300-6466 = {}".format(
                lab, r.chr1, r.id, r.length, len(inside), len(near)))
