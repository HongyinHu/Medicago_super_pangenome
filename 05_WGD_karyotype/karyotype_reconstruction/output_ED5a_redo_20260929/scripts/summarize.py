#!/usr/bin/env python3
"""Per-species filtering / mapping summary for the ED5a redo."""
import os

import pandas as pd

ROOT = "path/to/project/37.karyotype_reconstruction"
OUT = os.path.join(ROOT, "output_ED5a_redo_20260929")
W = os.path.join(OUT, os.environ.get("RUN", "02_wgdi"))
PREV = os.path.join(OUT, "02_wgdi")   # first filtered run (pvalue <= 0.05)

sp = pd.read_csv(os.path.join(OUT, "species_table.tsv"), sep="\t")
man = pd.read_csv(os.path.join(W, "run_manifest.tsv"), sep="\t")
sp = sp.merge(man[["label", "ks_cut"]], on="label")
rows = []
for _, s in sp.iterrows():
    d = os.path.join(W, s.label)
    raw = pd.read_csv(os.path.join(d, "in.blockinfo_raw"))
    c = pd.read_csv(os.path.join(d, s.label + "_aak.ortholog_blocks.csv"))
    o = pd.read_csv(os.path.join(d, s.label + "_aak.ortholog_blocks.1to1.csv"))
    lens = pd.read_csv(os.path.join(d, "in.lens1"), sep="\t", header=None, names=["chr", "bp", "genes"])
    km = pd.read_csv(os.path.join(d, "km_result.txt"), sep="\t", header=None,
                     names=["chr", "start", "end", "color", "cls"])
    old = pd.read_csv(os.path.join(ROOT, "output_Kary", s.subdir, "km_result.txt"), sep="\t", header=None,
                      names=["chr", "start", "end", "color", "cls"])
    prev = pd.read_csv(os.path.join(PREV, s.label, "km_result.txt"), sep="\t", header=None,
                       names=["chr", "start", "end", "color", "cls"])
    drop = pd.read_csv(os.path.join(d, s.label + "_aak.ortholog_blocks.1to1.csv.dropped.tsv"), sep="\t")
    total_genes = lens.genes.sum()
    colored = (km.end - km.start + 1).sum()
    prev_colored = (prev.end - prev.start + 1).sum()
    rows.append({
        "label": s.label, "subdir": s.subdir, "ks_cut": s.ks_cut,
        "chr_in_lens": len(lens), "chr_in_km": km.chr.nunique(),
        "blocks_raw": len(raw),
        "blocks_raw_ge20": int((raw.length >= 20).sum()),
        "blocks_after_c": len(c), "blocks_1to1": len(o),
        "genes_in_1to1_blocks": int(o.length.sum()),
        "blocks_dropped_1to1": len(drop),
        "segments_old": len(old), "segments_prev_run": len(prev), "segments_new": len(km),
        "uncolored_genes_prev_run": int(total_genes - prev_colored),
        "uncolored_genes_new": int(total_genes - colored),
        "frac_genome_genes_in_segments": round(colored / float(total_genes), 3),
        "amk_colors_used": km.color.nunique(),
        "status": "done" if os.path.exists(os.path.join(d, "status.done")) else "FAILED",
    })
df = pd.DataFrame(rows)
df.to_csv(os.path.join(W, "filter_summary.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(df.to_string(index=False))
