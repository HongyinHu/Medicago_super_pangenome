#!/usr/bin/env python3
"""Where are the AMK8 genes around CEN8 (AMK8 2350-2600) in M. praecox?
Best reciprocal-direction evidence: for each AMK8 gene, Mpra genes whose best hit it is."""
import pandas as pd

R = "path/to/project/37.karyotype_reconstruction"
b = pd.read_csv(R + "/output_Kary/genome_410/all.blastp.txt", sep="\t", header=None)
b = b[(b[10] <= 1e-5) & (b[11] >= 100)].sort_values([0, 11], ascending=[True, False]).drop_duplicates(0)
aak = pd.read_csv(R + "path/to/data", sep="\t", header=None, usecols=range(6),
                  names=["chr", "id", "s", "e", "str", "order"]).set_index("id")
pg = pd.read_csv(R + "/output_ED5a_redo_20260929/04_wgdi_p0.2/Mpra/in.gff1", sep="\t", header=None,
                 usecols=range(6), names=["chr", "gene", "s", "e", "str", "order"]).set_index("gene")
b = b[b[1].isin(aak.index) & b[0].isin(pg.index)]
b["achr"] = aak.loc[b[1], "chr"].values
b["aord"] = aak.loc[b[1], "order"].values
b["qchr"] = pg.loc[b[0], "chr"].values
b["qmb"] = pg.loc[b[0], "s"].values / 1e6
w = b[(b.achr == 8) & b.aord.between(2350, 2600)].sort_values("aord")
print("AMK8 genes 2350-2600 with an Mpra best hit: %d of 251" % w[1].nunique())
w["bin"] = (w.aord // 25) * 25
print(w.groupby(["bin", "qchr"]).agg(n=("qmb", "size"), mb_min=("qmb", "min"), mb_max=("qmb", "max")).to_string())
