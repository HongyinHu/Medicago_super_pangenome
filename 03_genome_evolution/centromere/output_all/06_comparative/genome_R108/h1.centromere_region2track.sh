python - <<'PY'
import pandas as pd

infile = "genome_R108.centromere_regions.tsv"
outfile = "genome_R108.chromosome_centromere_track.bed9"

df = pd.read_csv(infile, sep="\t")

# 颜色
grey = "180,180,180"          # 染色体臂
peri = "198,184,215"          # 周着丝粒，浅紫
cen = "128,105,155"           # CENH3着丝粒，深紫

records = []

for _, r in df.iterrows():
    chrom = r["chr"]
    size = int(r["chr_size"])
    ps = int(r["peri_start"])
    cs = int(r["cen_start"])
    ce = int(r["cen_end"])
    pe = int(r["peri_end"])
    label = str(r["label"])

    # 左臂
    if ps > 0:
        records.append([chrom, 0, ps, label, 0, ".", 0, ps, grey])

    # 左周着丝粒
    if cs > ps:
        records.append([chrom, ps, cs, label, 0, ".", ps, cs, peri])

    # CENH3功能着丝粒
    if ce > cs:
        records.append([chrom, cs, ce, label, 0, ".", cs, ce, cen])

    # 右周着丝粒
    if pe > ce:
        records.append([chrom, ce, pe, label, 0, ".", ce, pe, peri])

    # 右臂
    if size > pe:
        records.append([chrom, pe, size, label, 0, ".", pe, size, grey])

out = pd.DataFrame(records)
out.to_csv(outfile, sep="\t", header=False, index=False)

print(f"written: {outfile}")
PY