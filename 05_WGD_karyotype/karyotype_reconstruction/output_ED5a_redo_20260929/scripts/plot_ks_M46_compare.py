#!/usr/bin/env python3
"""Intra-genome Ks distribution of M. fischeriana: old (M46_1) vs new (M46_2) assembly.

Reads WGDI -kp outputs (<run>/Mfi_Mfi.kspeak.csv, i.e. the blocks passing the
[kspeaks] filters) and draws Gaussian KDEs of block-median, block-average and
all-pair Ks over 0-3, as in the WGDI kspeak plot. Also reports the WGD peak
positions so the two assemblies can be compared numerically.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

W = "path/to/project/5.WGD_2/output"
RUNS = [("old assembly (M46_1, 8 chr)", "X_genome_M46_1"), ("new assembly (M46_2, 7 chr)", "X_genome_M46_2")]
X = np.linspace(0, 3, 600)
matplotlib.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
                            "pdf.fonttype": 42})


def series(path):
    d = pd.read_csv(path)
    pairs = []
    for s in d["ks"].astype(str):
        pairs += [float(k) for k in s.split("_") if k not in ("", "nan")]
    pairs = np.array(pairs)
    med = d["ks_median"].astype(float).values
    avg = d["ks_average"].astype(float).values
    keep = lambda a: a[(a >= 0) & (a <= 3)]
    return {"block median": keep(med), "block average": keep(avg), "all pairs": keep(pairs)}, len(d)


def wgd_peak(x, y, lo=0.4, hi=1.5):
    m = (x >= lo) & (x <= hi)
    return x[m][np.argmax(y[m])]


fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
rows = []
for ax, (title, run) in zip(axes, RUNS):
    f = os.path.join(W, run, "Mfi_Mfi.kspeak.csv")
    if not os.path.exists(f):
        ax.set_title(title + "\n(not available)")
        continue
    S, nblk = series(f)
    for (lab, a), col in zip(S.items(), ["#E8272A", "black", "#2B5CB8"]):
        y = gaussian_kde(a)(X)
        ax.plot(X, y, color=col, lw=1.1, label=lab)
        rows.append(dict(assembly=run, series=lab, n=len(a), wgd_peak=round(wgd_peak(X, y), 3),
                         frac_le_0p3=round(float(np.mean(a <= 0.3)), 3)))
    ax.set_title("M. fischeriana  " + title, fontsize=10, style="italic")
    ax.set_xlim(0, 3)
    ax.set_xlabel("Ks", style="italic")
    ax.grid(alpha=0.4)
    ax.text(0.97, 0.55, "%d blocks" % nblk, transform=ax.transAxes, ha="right", fontsize=8)
axes[0].set_ylabel("Frequency")
axes[0].legend(fontsize=8, frameon=True)
fig.tight_layout()
out = os.path.join(W, "X_genome_M46_2")
fig.savefig(os.path.join(out, "Mfi_ks_old_vs_new.pdf"))
fig.savefig(os.path.join(out, "Mfi_ks_old_vs_new.png"), dpi=300)
r = pd.DataFrame(rows)
r.to_csv(os.path.join(out, "Mfi_ks_old_vs_new.summary.tsv"), sep="\t", index=False)
print(r.to_string(index=False))
