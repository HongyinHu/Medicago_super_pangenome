#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


base = Path("path/to/project/10.centromere_analysis/output_synthen/06_functional_centromere_sequence_composition/R108_TEsorter_TRASH_test")
tables = {
    "tandem_first": base / "genome_R108.functional_centromere_sequence_composition.broad.average.tsv",
    "te_first": base / "TE_first/genome_R108.functional_centromere_sequence_composition.broad.average.tsv",
}

frames = []
for mode, path in tables.items():
    df = pd.read_csv(path, sep="\t")
    df.insert(0, "priority_mode", mode)
    frames.append(df)

out = pd.concat(frames, ignore_index=True)
out.to_csv(base / "genome_R108.functional_centromere_sequence_composition.priority_mode_average_comparison.tsv", sep="\t", index=False, float_format="%.4f")
print(out.to_string(index=False))
