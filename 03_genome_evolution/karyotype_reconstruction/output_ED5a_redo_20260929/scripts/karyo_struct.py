#!/usr/bin/env python3
"""Per-genome chromosome composition in AMK terms, from the filtered
one-to-one ortholog blocks (04_wgdi_p0.2).

For each block: species chr/gene span, AMK chr/gene span, orientation
(+ if AMK order increases along the species chromosome).  Consecutive blocks
on the same species chromosome that continue the same AMK chromosome in the
same orientation are merged.  Output: 07_karyotype_evolution/segments.tsv and a
compact per-chromosome formula printed to stdout.
"""
import os
import sys

import numpy as np
import pandas as pd

OUT = ("path/to/project/37.karyotype_reconstruction/"
       "output_ED5a_redo_20260929")
W = os.path.join(OUT, "04_wgdi_p0.2")
DST = os.path.join(OUT, "07_karyotype_evolution")
MIN_GENES = int(sys.argv[1]) if len(sys.argv) > 1 else 30
MERGE_GAP = 400   # max AMK-coordinate jump (genes) still treated as the same segment

AMK_LEN = {1: 6320, 2: 5612, 3: 6726, 4: 6736, 5: 5544, 6: 5318, 7: 6466, 8: 6214}


def main():
    os.makedirs(DST, exist_ok=True)
    sp = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t").sort_values("order")
    rows = []
    for _, s in sp.iterrows():
        d = os.path.join(W, s.label)
        blk = pd.read_csv(os.path.join(d, s.label + "_aak.ortholog_blocks.1to1.csv"))
        lens = pd.read_csv(os.path.join(d, "in.lens1"), sep="\t", header=None,
                           names=["chr", "bp", "genes"])
        gff = pd.read_csv(os.path.join(d, "in.gff1"), sep="\t", header=None,
                          usecols=[0, 2, 3, 5], names=["chr", "s", "e", "order"])
        gff["chr"] = gff["chr"].astype(str)
        gff = gff.set_index(["chr", "order"])
        recs = []
        for _, r in blk.iterrows():
            b1 = np.array(list(map(int, r.block1.split("_"))))
            b2 = np.array(list(map(int, r.block2.split("_"))))
            if len(b1) < MIN_GENES:
                continue
            ori = "+" if np.corrcoef(b1, b2)[0, 1] >= 0 else "-"
            recs.append(dict(chr=str(r.chr1), s1=b1.min(), e1=b1.max(), amk=int(r.chr2),
                             s2=b2.min(), e2=b2.max(), ori=ori, n=len(b1)))
        df = pd.DataFrame(recs).sort_values(["chr", "s1"])
        merged = []
        for _, r in df.iterrows():
            r = r.to_dict()
            if merged:
                m = merged[-1]
                cont = (m["chr"] == r["chr"] and m["amk"] == r["amk"] and m["ori"] == r["ori"])
                if cont:
                    step = (r["s2"] - m["e2"]) if r["ori"] == "+" else (m["s2"] - r["e2"])
                    cont = -MERGE_GAP <= step <= MERGE_GAP
                if cont:
                    m["e1"] = max(m["e1"], r["e1"])
                    m["s2"], m["e2"] = min(m["s2"], r["s2"]), max(m["e2"], r["e2"])
                    m["n"] += r["n"]
                    continue
            merged.append(r)
        for m in merged:
            bp0 = gff.loc[(m["chr"], int(m["s1"])), "s"]
            bp1 = gff.loc[(m["chr"], int(m["e1"])), "e"]
            rows.append(dict(order=s.order, species=s.label, chr=m["chr"],
                             chr_genes=int(lens.set_index(lens.chr.astype(str)).loc[m["chr"], "genes"]),
                             start=int(m["s1"]), end=int(m["e1"]),
                             start_bp=int(bp0), end_bp=int(bp1), Mb=round((bp1 - bp0) / 1e6, 2),
                             AMK=m["amk"], amk_start=int(m["s2"]), amk_end=int(m["e2"]),
                             amk_frac=round((m["e2"] - m["s2"] + 1) / AMK_LEN[m["amk"]], 2),
                             ori=m["ori"], genes=int(m["n"])))
    seg = pd.DataFrame(rows)
    seg.to_csv(os.path.join(DST, "segments_min{}.tsv".format(MIN_GENES)), sep="\t", index=False)
    for (o, spc), g in seg.groupby(["order", "species"], sort=True):
        print("#", spc)
        for ch, c in g.groupby("chr", sort=False):
            parts = ["{}{}[{}-{}]{}".format("A", r.AMK, r.amk_start, r.amk_end, r.ori) +
                     "({}g,{}Mb)".format(r.genes, r.Mb) for _, r in c.iterrows()]
            print("  chr{:>2} ({}g): {}".format(ch, c.chr_genes.iloc[0], "  ".join(parts)))


if __name__ == "__main__":
    main()
