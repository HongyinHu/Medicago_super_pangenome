#!/usr/bin/env python
# fasta2phy.py
# Written by Nicolas D. Franco Sierra
# Universidad EAFIT, may/2015
# It converts an alignment in FASTA format to Phylip format
# Usage: fasta2phy.py <input_alignment.fasta> <output_alignment.phy>

import sys
from Bio import SeqIO

def length_counter(fasta):
    for seq in SeqIO.parse(fasta,"fasta"):
        return len(str(seq.seq))
        break

def taxa_counter(fasta):
    num_taxa = 0
    for record in SeqIO.parse(open(fasta, "rU"), "fasta"):
        num_taxa += 1
    return num_taxa


input_file = sys.argv[1]
output_file = sys.argv[2]

o = open(output_file,"w")


n_taxa = taxa_counter(input_file)
matrix_len = length_counter(input_file)

o.write(str(n_taxa) + " " + str(matrix_len) + "\n")

for record in SeqIO.parse(open(input_file, "rU"), "fasta"):
    o.write(str(record.id) + "   " + str(record.seq) + "\n")