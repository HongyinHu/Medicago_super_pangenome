import pandas as pd
import sys

sp = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

df = pd.read_csv(infile)
df = df.loc[:, ~df.columns.str.contains("^Unnamed|^$")]

records = []

for i, row in df.iterrows():
    chrom = str(row["name"])
    start = int(row["start"])
    end = int(row["end"]) + 1

    monomer_len = int(row["most.freq.value.N"])
    array_len = end - start

    rep_class = str(row["class"]) if "class" in df.columns else "NA"
    if rep_class == "NA" or rep_class == "nan" or rep_class == "":
        repeat_family = f"{sp}_TRASH_{monomer_len}bp_{i+1:06d}"
    else:
        repeat_family = rep_class

    records.append([
        chrom,
        start,
        end,
        repeat_family,
        monomer_len,
        array_len
    ])

pd.DataFrame(records).to_csv(outfile, sep="\t", header=False, index=False)
