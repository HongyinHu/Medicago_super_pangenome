#!/usr/bin/env python3
"""Inter-AMK junctions (adjacent major segments from different AMK chromosomes)
in every genome, and whether any junction is shared between genomes
(same AMK pair, both AMK breakpoints within TOL genes).  Also reports, for the
x = 7 genomes, which chromosomes carry more than one AMK and their size rank."""
import os
from itertools import combinations

import pandas as pd

OUT = ("path/to/project/37.karyotype_reconstruction/"
       "output_ED5a_redo_20260929")
D = os.path.join(OUT, "07_karyotype_evolution")
MIN_GENES, TOL = 100, 300
X7 = ["Mpra", "Mpol", "Mfis"]

seg = pd.read_csv(os.path.join(D, "segments_min30.tsv"), sep="\t")
seg["chr"] = seg["chr"].astype(str)
seg = seg[seg.genes >= MIN_GENES]
J = []
for (sp, ch), g in seg.groupby(["species", "chr"]):
    g = g.sort_values("start").to_dict("records")
    # merge consecutive pieces of the same AMK
    m = []
    for r in g:
        if m and m[-1]["AMK"] == r["AMK"]:
            m[-1]["end"] = r["end"]
            m[-1]["last"] = r
        else:
            r["first"], r["last"] = r, r
            m.append(r)
    for a, b in zip(m, m[1:]):
        la, fb = a["last"], b["first"]
        bpa = la["amk_end"] if la["ori"] == "+" else la["amk_start"]
        bpb = fb["amk_start"] if fb["ori"] == "+" else fb["amk_end"]
        J.append(dict(species=sp, chr=ch, amk_a=la["AMK"], bp_a=bpa, amk_b=fb["AMK"], bp_b=bpb,
                      species_pos=int(la["end"])))
J = pd.DataFrame(J)
J.to_csv(os.path.join(D, "inter_AMK_junctions.tsv"), sep="\t", index=False)
print("junctions per genome:", J.groupby("species").size().to_dict())


def key(r):
    return (min(r.amk_a, r.amk_b), max(r.amk_a, r.amk_b))


def bps(r):
    return {r.amk_a: r.bp_a, r.amk_b: r.bp_b}


shared = []
for (i, r), (j, s) in combinations(J.iterrows(), 2):
    if r.species == s.species or key(r) != key(s):
        continue
    br, bs = bps(r), bps(s)
    if all(abs(br[k] - bs[k]) <= TOL for k in br):
        shared.append((r.species, r.chr, s.species, s.chr, key(r), br, bs))
print("\nshared junctions (same AMK pair, both breakpoints within {} genes):".format(TOL))
for x in shared:
    print("  {}:chr{} ~ {}:chr{}  AMK{}  {} vs {}".format(*x))
if not shared:
    print("  none")

print("\njunctions of the x = 7 genomes and the closest junction of the same AMK pair elsewhere:")
for sp in X7:
    for _, r in J[J.species == sp].iterrows():
        oth = J[(J.species != sp)]
        oth = oth[oth.apply(lambda s: key(s) == key(r), axis=1)]
        best = ""
        if len(oth):
            d = oth.apply(lambda s: max(abs(bps(s)[k] - bps(r)[k]) for k in bps(r)), axis=1)
            k = d.idxmin()
            best = "{}:chr{} (max bp diff {} genes)".format(oth.loc[k, "species"], oth.loc[k, "chr"], int(d[k]))
        print("  {} chr{}: AMK{}@{} | AMK{}@{}   nearest same-pair: {}".format(
            sp, r.chr, r.amk_a, r.bp_a, r.amk_b, r.bp_b, best or "-"))

print("\nchromosome size ranks (genes) of the x = 7 genomes; * = carries >1 AMK")
for sp in X7:
    g = seg[seg.species == sp].groupby("chr").agg(genes=("chr_genes", "first"), amks=("AMK", "nunique"))
    g = g.sort_values("genes", ascending=False)
    print("  " + sp + ": " + ", ".join("chr{}{}={}".format(c, "*" if r.amks > 1 else "", r.genes)
                                     for c, r in g.iterrows()))
