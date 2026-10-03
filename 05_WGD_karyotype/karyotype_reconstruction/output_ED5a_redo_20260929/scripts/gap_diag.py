#!/usr/bin/env python3
"""List every internal gap between km_result segments and explain each with the
raw WGDI blocks that overlap it and the filter step that removed them."""
import os
import sys

import pandas as pd

W = ("path/to/project/37.karyotype_reconstruction/"
     "output_ED5a_redo_20260929/" + os.environ.get("RUN", "02_wgdi"))
MANIFEST = pd.read_csv(os.path.join(W, "run_manifest.tsv"), sep="\t").set_index("label")


def span(r):
    b = [int(x) for x in str(r["block1"]).split("_") if x]
    return min(b), max(b)


def main(label, min_gap):
    d = os.path.join(W, label)
    cut = MANIFEST.loc[label, "ks_cut"]
    pcut = MANIFEST.loc[label, "pvalue"] if "pvalue" in MANIFEST.columns else 0.05
    km = pd.read_csv(os.path.join(d, "km_result.txt"), sep="\t", header=None,
                     names=["chr", "start", "end", "color", "cls"])
    raw = pd.read_csv(os.path.join(d, "in.blockinfo_raw"))
    c = set(pd.read_csv(os.path.join(d, label + "_aak.ortholog_blocks.csv"))["id"])
    o = set(pd.read_csv(os.path.join(d, label + "_aak.ortholog_blocks.1to1.csv"))["id"])
    raw[["s1", "e1"]] = raw.apply(lambda r: pd.Series(span(r)), axis=1)
    for col in ("chr1", "chr2"):
        raw[col] = raw[col].astype(str)
    km["chr"] = km["chr"].astype(str)
    print("#### {}  (Ks cut {})".format(label, cut))
    for ch, g in km.groupby("chr"):
        g = g.sort_values("start").reset_index(drop=True)
        for i in range(1, len(g)):
            gs, ge = g.loc[i - 1, "end"] + 1, g.loc[i, "start"] - 1
            if ge - gs + 1 < min_gap:
                continue
            print("\n== chr{} genes {}-{} ({} genes) between {} and {}".format(
                ch, gs, ge, ge - gs + 1, g.loc[i - 1, "color"], g.loc[i, "color"]))
            ov = raw[(raw.chr1 == ch) & (raw.s1 <= ge) & (raw.e1 >= gs)].copy()
            if ov.empty:
                print("   no collinear block with AMK at all (-icl found nothing)")
                continue
            ks_frac = []
            for _, r in ov.iterrows():
                ks = [float(k) for k in str(r["ks"]).split("_") if k not in ("", "nan")]
                ks_frac.append(sum(0 <= k <= cut for k in ks) / float(len(ks)) if ks else 0)
            ov["frac_in_peak"] = ks_frac

            def why(r):
                if r["id"] in o:
                    return "KEPT (colour run <20 genes or overridden)"
                if r["id"] in c:
                    return "dropped: one-to-one"
                if r["length"] < 20:
                    return "dropped: <20 genes"
                if r["pvalue"] > pcut:
                    return "dropped: pvalue"
                if r["frac_in_peak"] < 0.5:
                    return "dropped: Ks (WGD paralog)"
                return "dropped: -c overlap/homo"
            ov["fate"] = ov.apply(why, axis=1)
            ov["ovl_genes"] = ov.apply(lambda r: min(r.e1, ge) - max(r.s1, gs) + 1, axis=1)
            print(ov.sort_values("length", ascending=False)[
                ["id", "chr2", "s1", "e1", "length", "ovl_genes", "ks_median",
                 "frac_in_peak", "pvalue", "fate"]].head(12).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
