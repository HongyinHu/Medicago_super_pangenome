# 先获取泛基因组类别的基因列表，后根据基因列表获取序列

import re
import os
import sys
from Bio import SeqIO

def main(genelist_file,all_pep_fa):
	genefam_list = list()
	with open(genelist_file) as inf:
		for lines in inf.readlines():
			genefam_list.append(lines.strip())

	gene_pep = dict()
	for seq in SeqIO.parse(all_pep_fa,"fasta"):
		gene_pep[seq.id]=str(seq.seq)

	for i in genefam_list:
		print(f">{i}\n{gene_pep[i]}")

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s genelist_file all_pep_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])



