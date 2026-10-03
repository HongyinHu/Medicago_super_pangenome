# 先获取泛基因组类别的基因列表，后根据基因列表获取序列

import re
import os
import sys
from Bio import SeqIO

def main(genefam_file,Orthogroups_tsv,all_pep_fa):
	genefam_list = list()
	with open(genefam_file) as inf:
		for lines in inf.readlines():
			genefam_list.append(lines.strip())


	genefam_gene = dict()
	with open(Orthogroups_tsv) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index > 0:
				line = lines.strip().split("\t")
				genelist = list()
				for i in line[1:]:
					gene_id_list = re.split(r",",i)
					genelist += gene_id_list
				genefam_gene[line[0]] = genelist

	gene_pep = dict()
	for seq in SeqIO.parse(all_pep_fa,"fasta"):
		gene_pep[seq.id]=str(seq.seq)

	with open(f"{os.path.basename(genefam_file)}.genelist","w") as ouf1:
		for i in genefam_list:
			if i in genefam_gene:
				for ii in genefam_gene[i]:
					ii = ii.strip()
					if ii in gene_pep:
						print(ii,file=ouf1)

	with open(f"{os.path.basename(genefam_file)}.genelist.pep","w") as ouf2:
		for i in genefam_list:
			if i in genefam_gene:
				for ii in genefam_gene[i]:
					ii = ii.strip()
					if ii in gene_pep:
						print(f">{ii}\n{gene_pep[ii]}",file=ouf2)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s genefam_file Orthogroups_tsv all_pep_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3])



