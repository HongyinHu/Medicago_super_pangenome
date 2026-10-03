#根据cafe5的输出结果，筛选某一个物种的快速扩张或快速收缩的基因家族
import re
import os
import sys
import pandas as pd

# ========= 需要改的参数 =========
species = "Mla_454<8>"   # 改成你的目标物种列名
species_out = "Mla_454"
alpha = 0.05             # 你想用的阈值，可改成 0.01
use_family_sig = True    # 是否再叠加 family-level 的 y 过滤
# ==============================

def read_table(path, na_values=None):
	try:
		return pd.read_csv(path, sep="\t", na_values=na_values)
	except Exception:
		return pd.read_csv(path, seq=r"\s+", engine="python", na_values=na_values)

fam = read_table("Gamma_family_results.txt")
chg = read_table("Gamma_change.tab")
bp  = read_table("Gamma_branch_probabilities.tab", na_values=["N/A"])

# family_results 第一列一般是 FamilyID / #FamilyID
fam_id_col = fam.columns[0]
chg_id_col = chg.columns[0]
bp_id_col  = bp.columns[0]

# 找到 family-level 显著列
sig_col_candidates = [c for c in fam.columns if "Significant" in str(c)]
if len(sig_col_candidates) == 0:
    raise ValueError("在 Gamma_family_results.txt 中没有找到 'Significant' 列。")
sig_col = sig_col_candidates[0]

# 统一列名
fam = fam.rename(columns={fam_id_col: "FamilyID", sig_col: "family_sig"})
chg = chg.rename(columns={chg_id_col: "FamilyID", species: "change"})
bp  = bp.rename(columns={bp_id_col: "FamilyID", species: "prob"})

# 只保留要用的列
chg = chg[["FamilyID", "change"]].copy()
bp  = bp[["FamilyID", "prob"]].copy()
fam = fam[["FamilyID", "family_sig"]].copy()

# 转数值
chg["change"] = pd.to_numeric(chg["change"], errors="coerce")
bp["prob"] = pd.to_numeric(bp["prob"], errors="coerce")

# 转换为字符串
chg["FamilyID"] = chg["FamilyID"].astype(str)
bp["FamilyID"] = bp["FamilyID"].astype(str)
fam["FamilyID"] = fam["FamilyID"].astype(str)

# 合并
df = chg.merge(bp, on="FamilyID", how="inner").merge(fam, on="FamilyID", how="left")

# family-level 过滤（可选）
if use_family_sig:
    df = df[df["family_sig"].astype(str).str.lower().eq("y")]

# species-specific 提取
expanded = df[(df["prob"] < alpha) & (df["change"] > 0)].copy()
contracted = df[(df["prob"] < alpha) & (df["change"] < 0)].copy()

# 输出 family ID 列表
expanded[["FamilyID"]].to_csv(f"{species_out}_expanded_families.txt", index=False, header=False)
contracted[["FamilyID"]].to_csv(f"{species_out}_contracted_families.txt", index=False, header=False)

# 同时输出带证据的表，方便你核对
expanded.to_csv(f"{species_out}_expanded_families.full.tsv", sep="\t", index=False)
contracted.to_csv(f"{species_out}_contracted_families.full.tsv", sep="\t", index=False)

print(f"[{species}] 显著扩张 family 数量:", len(expanded))
print(f"[{species}] 显著收缩 family 数量:", len(contracted))
print("已输出：")
print(f"  {species}_expanded_families.txt")
print(f"  {species}_contracted_families.txt")
print(f"  {species}_expanded_families.full.tsv")
print(f"  {species}_contracted_families.full.tsv")