
from collections import Counter

aln = "candidate_array.monomers.sample200.aln.fa"
out = "candidate_array.final_monomer_consensus.fa"

seqs = []
name = None
seq = []

with open(aln) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs.append("".join(seq).upper())
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    if name is not None:
        seqs.append("".join(seq).upper())

length = max(len(s) for s in seqs)
consensus = []

for i in range(length):
    bases = []
    for s in seqs:
        if i < len(s):
            b = s[i]
            if b in "ACGT":
                bases.append(b)
    if not bases:
        continue
    c = Counter(bases)
    base, count = c.most_common(1)[0]
    if count / len(bases) >= 0.5:
        consensus.append(base)
    else:
        consensus.append("N")

consensus = "".join(consensus)

with open(out, "w") as f:
    f.write(">final_monomer_consensus\n")
    for i in range(0, len(consensus), 80):
        f.write(consensus[i:i+80] + "\n")

print(f"consensus length: {len(consensus)}")
print(f"written: {out}")