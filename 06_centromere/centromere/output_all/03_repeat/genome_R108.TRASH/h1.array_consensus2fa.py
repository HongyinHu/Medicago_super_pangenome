import pandas as pd


SP="genome_R108"

infile = f"Summary.of.repetitive.regions.{SP}.fa.csv"
outfile = f"{SP}.TRASH_candidate_consensus.fa"

df = pd.read_csv(infile)
df = df.loc[:, ~df.columns.str.contains("^Unnamed|^$")]

with open(outfile, "w") as out:
    for i, row in df.iterrows():
        chrom = str(row["name"])
        start = int(row["start"])
        end = int(row["end"])
        monomer_len = int(row["most.freq.value.N"])
        seq = str(row["consensus.primary"]).replace(" ", "").replace("-", "").upper()

        if seq == "NA" or seq == "NAN" or len(seq) < 20:
            continue

        name = f"TRASH_{i+1:06d}|{chrom}:{start}-{end}|len{monomer_len}|array{end-start+1}"
        out.write(f">{name}\n{seq}\n")

print(f"written: {outfile}")