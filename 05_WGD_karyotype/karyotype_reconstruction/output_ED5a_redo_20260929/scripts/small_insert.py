#!/usr/bin/env python3
"""Find short foreign-AMK segments embedded inside a chromosome (flanked on both
sides by the same other AMK) in every genome, with the WGDI blocks behind them."""
import os

import pandas as pd

OUT = ("path/to/project/37.karyotype_reconstruction/"
       "output_ED5a_redo_20260929")
W = os.path.join(OUT, os.environ.get("RUN", "04_wgdi_p0.2"))
AMK = {"royalblue": 1, "red": 2, "#99cc00": 3, "deepskyblue": 4, "#339966": 5,
       "#ffcc00": 6, "fuchsia": 7, "#aa6e20": 8}

sp = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t").sort_values("order")
rows = []
for _, s in sp.iterrows():
    d = os.path.join(W, s.label)
    km = pd.read_csv(os.path.join(d, "km_result.txt"), sep="\t", header=None,
                     names=["chr", "start", "end", "color", "cls"])
    km["amk"] = [AMK[c.lower()] for c in km.color]
    gff = pd.read_csv(os.path.join(d, "in.gff1"), sep="\t", header=None, usecols=[0, 1, 2, 3, 5],
                      names=["chr", "gene", "s", "e", "order"])
    gff["chr"] = gff["chr"].astype(str)
    gff = gff.set_index(["chr", "order"])
    blk = pd.read_csv(os.path.join(d, s.label + "_aak.ortholog_blocks.1to1.csv"))
    for col in ("chr1", "chr2"):
        blk[col] = blk[col].astype(str)
    for ch, g in km.groupby("chr"):
        g = g.sort_values("start").reset_index(drop=True)
        for i in range(1, len(g) - 1):
            a, m, b = g.loc[i - 1], g.loc[i], g.loc[i + 1]
            if a.amk == b.amk and m.amk != a.amk:
                bp0 = gff.loc[(str(ch), int(m.start)), "s"]
                bp1 = gff.loc[(str(ch), int(m.end)), "e"]
                # blocks that lie inside the inserted segment
                sub = blk[(blk.chr1 == str(ch)) & (blk.chr2 == str(m.amk))]
                sub = sub[sub.apply(lambda r: m.start <= min(map(int, r.block1.split("_"))) <= m.end, axis=1)]
                rows.append({
                    "species": s.label, "chr": ch, "host_AMK": int(a.amk), "insert_AMK": int(m.amk),
                    "gene_start": int(m.start), "gene_end": int(m.end),
                    "genes": int(m.end - m.start + 1), "Mb": round((bp1 - bp0) / 1e6, 2),
                    "blocks": ";".join("id{}:{}g,Ks{:.3f},AMK{}:{}-{}".format(
                        r.id, r.length, r.ks_median, r.chr2,
                        min(map(int, r.block2.split("_"))), max(map(int, r.block2.split("_"))))
                        for _, r in sub.iterrows()),
                })
df = pd.DataFrame(rows)
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 120)
print(df.to_string(index=False))
df.to_csv(os.path.join(OUT, "05_figure_p0.2", "embedded_small_segments.tsv"), sep="\t", index=False)
