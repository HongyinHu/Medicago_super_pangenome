#!/usr/bin/env python3
"""Survey block-level and pair-level Ks of every species-vs-AMK WGDI run.

Reads the existing WGDI [blockinfo] outputs (``-icl`` -> ``-ks`` -> ``-bi``)
and writes per-block tables plus a multi-panel Ks histogram so that the
ortholog (speciation) peak and the Papilionoideae WGD peak can be separated.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "path/to/project/37.karyotype_reconstruction"

# (short label, full species name, output_Kary subdirectory)
SPECIES = [
    ("Malbus", "Melilotus albus", "genome_Mal"),
    ("Marc", "Medicago arborea", "genome_Mar"),
    ("Mrut", "Medicago ruthenica", "genome_Mru"),
    ("Mlan", "Medicago lanigera", "genome_454"),
    ("Mcar", "Medicago carstiensis", "genome_474"),
    ("Mcre", "Medicago cretacea", "genome_468"),
    ("Msat_cae", "Medicago sativa subsp. caerulea", "genome_Msa"),
    ("Msat_zm4", "Medicago sativa ZM4", "genome_ZM4"),
    ("Mmar", "Medicago marina", "genome_457"),
    ("Mpra", "Medicago praecox", "genome_410"),
    ("Mtru_R108", "Medicago truncatula R108", "genome_R108"),
    ("Mpol", "Medicago polymorpha", "genome_Mpo_rerun_20260618"),
    ("Morb", "Medicago orbicularis", "genome_M22"),
    ("Msec", "Medicago secundiflora", "genome_461"),
    ("Mlup", "Medicago lupulina", "genome_395"),
    ("Msuf", "Medicago suffruticosa", "genome_472"),
    ("Medg", "Medicago edgeworthii", "genome_482"),
    ("Mfis", "Medicago fischeriana", "genome_M46_2"),
    ("Mrad", "Medicago radiata", "genome_436"),
]


def parse_section(path, wanted):
    sec, out = None, {}
    for raw in open(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            sec = line[1:-1].strip()
            continue
        if sec == wanted and "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def split_ks(text):
    vals = []
    for k in str(text).split("_"):
        if k == "" or k == "nan":
            continue
        try:
            vals.append(float(k))
        except ValueError:
            pass
    return vals


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    block_rows, peak_rows = [], []
    fig, axes = plt.subplots(5, 4, figsize=(20, 20))
    axes = axes.ravel()
    bins = np.arange(0, 2.51, 0.02)
    for i, (label, full, sub) in enumerate(SPECIES):
        d = os.path.join(ROOT, "output_Kary", sub)
        bi = parse_section(os.path.join(d, "total.conf"), "blockinfo")
        path = os.path.join(d, bi["savefile"])
        df = pd.read_csv(path)
        pair_ks = []
        for _, r in df.iterrows():
            ks = split_ks(r["ks"])
            ks = [k for k in ks if k >= 0]
            pair_ks.extend(ks)
            block_rows.append({
                "species": label, "block_id": r["id"], "chr1": r["chr1"], "chr2": r["chr2"],
                "length": r["length"], "pvalue": r["pvalue"], "homo1": r["homo1"],
                "ks_median": r["ks_median"], "n_ks": len(ks),
                "frac_ks_le_0.3": np.mean([k <= 0.3 for k in ks]) if ks else np.nan,
                "frac_ks_le_0.4": np.mean([k <= 0.4 for k in ks]) if ks else np.nan,
                "frac_ks_le_0.5": np.mean([k <= 0.5 for k in ks]) if ks else np.nan,
            })
        pk = np.array(pair_ks)
        h, edges = np.histogram(pk, bins=bins)
        centers = (edges[:-1] + edges[1:]) / 2
        # smooth with a 5-bin moving average before locating peaks/trough
        sm = np.convolve(h, np.ones(5) / 5, mode="same")
        lo = centers <= 0.45
        hi = (centers > 0.45) & (centers <= 1.3)
        p1 = centers[lo][np.argmax(sm[lo])]
        p2 = centers[hi][np.argmax(sm[hi])]
        between = (centers > p1) & (centers < p2)
        trough = centers[between][np.argmin(sm[between])] if between.any() else np.nan
        peak_rows.append({"species": label, "full_name": full, "subdir": sub,
                          "n_blocks": len(df), "n_pairs": len(pk),
                          "ortholog_peak_ks": round(p1, 3), "wgd_peak_ks": round(p2, 3),
                          "trough_ks": round(trough, 3),
                          "frac_pairs_le_trough": round(float(np.mean(pk <= trough)), 3)})
        ax = axes[i]
        ax.bar(centers, h, width=0.02, color="#5B7DB1")
        ax.plot(centers, sm, color="k", lw=0.8)
        ax.axvline(trough, color="red", ls="--", lw=1)
        ax.set_yscale("log")
        ax.set_xlim(0, 2.5)
        ax.set_title("{} ({})  peak={:.2f} WGD={:.2f} cut={:.2f}".format(label, sub, p1, p2, trough), fontsize=9)
    for j in range(len(SPECIES), len(axes)):
        axes[j].axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "ks_pair_histograms_log.png"), dpi=80)
    pd.DataFrame(block_rows).to_csv(os.path.join(outdir, "blocks_ks_survey.tsv"), sep="\t", index=False)
    pd.DataFrame(peak_rows).to_csv(os.path.join(outdir, "ks_peaks_by_species.tsv"), sep="\t", index=False)
    print(pd.DataFrame(peak_rows).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
