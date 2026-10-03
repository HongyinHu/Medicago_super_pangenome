#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

out = Path("path/to/project/N_4.pod_spiny/05_direction_free_sv_cosegregation_20260701")

a = pd.read_csv(out / "results/direction_free_all_candidates.tsv", sep="\t")
b = pd.read_csv(out / "results/direction_free_igv_selected.tsv", sep="\t")
c = pd.read_csv(out / "results/candidate_strict_flankQC_summary.tsv", sep="\t")

print("ALL by ref/SVTYPE")
print("columns:", ",".join(a.columns[:40]))
ref_col = "ref" if "ref" in a.columns else "igv_ref" if "igv_ref" in a.columns else "best_ref" if "best_ref" in a.columns else "dataset"
print(a.groupby([ref_col, "SVTYPE"]).size().to_string())
print()
print("IGV selected by ref/SVTYPE/best_direction")
b_ref_col = "ref" if "ref" in b.columns else "igv_ref"
print(b.groupby([b_ref_col, "SVTYPE", "best_direction"]).size().to_string())
print()
print("flankQC by ref/SVTYPE/verdict/direction")
print(c.groupby(["ref", "SVTYPE", "strict_flankQC_verdict", "association_direction"]).size().to_string())
