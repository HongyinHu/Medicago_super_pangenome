#!/usr/bin/env python3
import sys
from pathlib import Path

if len(sys.argv) not in (2, 3):
    raise SystemExit("Usage: pad_fasta_alignment.py in.fa [out.fa]")
infile = Path(sys.argv[1])
outfile = Path(sys.argv[2]) if len(sys.argv) == 3 else infile
records = []
name = None
seq_parts = []
for line in infile.read_text().splitlines():
    if not line:
        continue
    if line.startswith(">"):
        if name is not None:
            records.append((name, "".join(seq_parts)))
        name = line[1:].strip().split()[0]
        seq_parts = []
    else:
        seq_parts.append(line.strip())
if name is not None:
    records.append((name, "".join(seq_parts)))
if not records:
    raise SystemExit(f"No FASTA records in {infile}")
max_len = max(len(seq) for _, seq in records)
with outfile.open("w") as out:
    for name, seq in records:
        padded = seq + ("-" * (max_len - len(seq)))
        out.write(f">{name}\n")
        for i in range(0, len(padded), 80):
            out.write(padded[i:i+80] + "\n")
print(f"padded {len(records)} records to {max_len} columns: {outfile}")
