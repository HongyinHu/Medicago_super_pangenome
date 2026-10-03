#!/usr/bin/env python3
import numpy as np
import pandas as pd


prefixes = [
    "functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.similarity_matrix.tsv",
    "functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.similarity_matrix.tsv",
    "functional_centromere.intact_LTRRT.similarity_heatmap.Copia_GypsyTop5.BLASTN_maxPident.minCov20.similarity_matrix.tsv",
]

for prefix in prefixes:
    matrix = pd.read_csv(prefix, sep="\t", index_col=0).to_numpy(float)
    n = matrix.shape[0]
    values = matrix[np.triu_indices(n, 1)]
    nonzero = values[values > 0]
    print(prefix)
    print(f"  pairs: {values.size}")
    print(f"  nonzero: {nonzero.size}")
    print(f"  nonzero_pct: {nonzero.size / values.size * 100:.2f}")
    if nonzero.size:
        q = np.percentile(nonzero, [50, 75, 90, 95, 100])
        print("  p50/p75/p90/p95/max: " + ", ".join(f"{x:.3f}" for x in q))
