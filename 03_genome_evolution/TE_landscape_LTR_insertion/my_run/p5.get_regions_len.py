#获取区间长度累加之和
import re
import os
import sys
from Bio import SeqIO

def main(file,genome_fa):
	gene_len = 0
	for seq in SeqIO.parse(genome_fa,"fasta"):
		gene_len += len(seq.seq)

	seq_len = 0
	with open(file) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			lens = int(line[2])-int(line[1])
			seq_len += lens

	print(genome_fa,seq_len,seq_len/gene_len)

if __name__ == '__main__':
	if len(sys.argv) !=  3:
		print("python %s bed_file genome_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])
