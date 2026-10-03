# 获取每个物种在泛基因组中的各种类型的基因表达水平

import re
import os
import sys

def species_gene_type(file,species_prefix,gene_tpm,type):
	gene_list = list()
	with open(file) as inf, open(f"{os.path.basename(file)}.{species_prefix}.{type}.tpm","w") as ouf:
		for lines in inf.readlines():
			line = lines.strip()
			if re.search(species_prefix,line):
				if line in gene_tpm:
					print("\t".join([type,gene_tpm[line]]),file=ouf)

def main(core_file,softcore_file,shell_file,specific_file,species_tpm,species_prefix):
	gene_tpm = dict()
	with open(species_tpm) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			gene_tpm["|".join([species_prefix,line[0]])] = line[1]

	species_gene_type(core_file,species_prefix,gene_tpm,"core")
	species_gene_type(softcore_file,species_prefix,gene_tpm,"softcore")
	species_gene_type(shell_file,species_prefix,gene_tpm,"shell")
	species_gene_type(specific_file,species_prefix,gene_tpm,"specific")


if __name__ == '__main__':
	if len(sys.argv) != 7:
		print("python %s core_file softcore_file shell_file specific_file species_tpm species_prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])

