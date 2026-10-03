#从获得基因组家族中查找具体的geneID

import re
import os
import sys
import pandas as pd

# ======== 需要改的参数 ========
species_family = "Mla_454<8>"   # CAFE里的列名
species_og = "Mla_454"          # Orthogroups.tsv里的物种列名，通常没有 <1>
# =============================

og = pd.read_csv("Orthogroups.tsv", sep="\t")

exp_ids = set(pd.read_csv(f"{species_og}_expanded_families.txt", header=None)[0].astype(str))
con_ids = set(pd.read_csv(f"{species_og}_contracted_families.txt", header=None)[0].astype(str))

family_col = og.columns[0]  # 一般是 Orthogroup

def split_genes(cell):
    if pd.isna(cell):
        return []
    s = str(cell).strip()
    if s == "":
        return []
    return [x.strip() for x in s.split(",") if x.strip()]

expanded_genes = []
for _, row in og[og[family_col].isin(exp_ids)].iterrows():
    expanded_genes.extend(split_genes(row[species_og]))

contracted_genes = []
for _, row in og[og[family_col].isin(con_ids)].iterrows():
    contracted_genes.extend(split_genes(row[species_og]))

with open(f"{species_og}_expanded_genes.txt", "w") as f:
    for g in expanded_genes:
        f.write(g + "\n")

with open(f"{species_og}_contracted_genes.txt", "w") as f:
    for g in contracted_genes:
        f.write(g + "\n")

print("expanded genes:", len(expanded_genes))
print("contracted genes:", len(contracted_genes))