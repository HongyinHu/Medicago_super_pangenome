#获取全长完整的LTR从genome.fa.mod.LTR.intact.gff3文件中

import re
import os
import sys
from Bio import SeqIO

def read_gff(file):
	gff_list = list()
	with open(file) as inf:
		for lines in inf.readlines():
			if not lines.isspace():
				line = lines.strip().split()
				name = re.search(r"ID=(.*?);",lines).group(1)
				gff_list.append([line[0],line[3],line[4],name,line[6]])

	return gff_list

def main(intact_gff,prefix):
	types = dict()

	with open(intact_gff) as inf, open(f"{prefix}.lLTR.gff","w") as ouf1, open(f"{prefix}.rLTR.gff","w") as ouf2:
		for lines in inf.readlines():
			if re.search("long_terminal_repeat",lines):
				if re.search("lLTR",lines):
					print(lines,file=ouf1)
				if re.search("rLTR",lines):
					print(lines,file=ouf2)

				clss = re.search(r"Classification=(.*?);",lines).group(1)
				types.setdefault(clss,[]).append(clss)

	lLTR_list = read_gff(f"{prefix}.lLTR.gff")
	rLTR_list = read_gff(f"{prefix}.rLTR.gff")

	with open(f"{prefix}.lLTR.bed","w") as ouf3:
		for i in lLTR_list:
			print("\t".join(i),file=ouf3)

	with open(f"{prefix}.rLTR.bed","w") as ouf4:
		for i in rLTR_list:
			print("\t".join(i),file=ouf4)

	for k,v in types.items():
		print(k,len(v))

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s intact_gff prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])







