import sys
from collections import Counter

aln = sys.argv[1]
out = sys.argv[2]
name = sys.argv[3]

seqs = []
seq = []

with open(aln) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq:
                seqs.append("".join(seq).upper())
            seq = []
        else:
            seq.append(line)
    if seq:
        seqs.append("".join(seq).upper())

if not seqs:
    raise SystemExit("No sequences found in alignment")

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

consensus = "".join(consensus).replace("-", "")

with open(out, "w") as f:
    f.write(f">{name}|len={len(consensus)}|ncopy={len(seqs)}\n")
    for i in range(0, len(consensus), 80):
        f.write(consensus[i:i+80] + "\n")

print(f"{name}\tconsensus_len={len(consensus)}\tncopy={len(seqs)}")
