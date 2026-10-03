# 对同源基因家族进行分类

import re
import os 
import sys
from Bio import SeqIO

def main(GeneCount,pep_fa,prefix):
	total_genes_list = set()

	for seq in SeqIO.parse(pep_fa,"fasta"):
		total_genes_list.add(seq.id)

	speclist_list = list()
	Orthogroup_list = list()
	Orthogroup_cont = dict()
	Orthogroup_dict = dict()

	with open(GeneCount) as inf:
		for index, lines in enumerate(inf.readlines()):
			line = lines.strip().split()
			if index == 0:
				speclist_list += line[1:]
			else:
				Orthogroup_list.append(line[0])
				Orthogroup_cont[line[0]] = list(map(int,line[1:-1]))

				for inde, i in enumerate(line[1:]):
					if line[0] not in Orthogroup_dict:
						Orthogroup_dict[line[0]] = {}
					Orthogroup_dict[line[0]][speclist_list[inde]] = int(i)

	single_copy_gene = 0
	multiple_copy_gene = 0
	unique_paralogs = 0
	Other_orthology_gene = 0

	do_Orthology = list()

	if prefix in speclist_list:
		# single_copy orthology
		for i in Orthogroup_list:
			if Orthogroup_dict[i][prefix] == 1 and all(x == 1 for x in Orthogroup_cont[i]):
				single_copy_gene += Orthogroup_dict[i][prefix]
				do_Orthology.append(i)

			if Orthogroup_dict[i][prefix] > 1 and all(x > 1 for x in Orthogroup_cont[i]):
				multiple_copy_gene += Orthogroup_dict[i][prefix]
				do_Orthology.append(i)

			if Orthogroup_dict[i][prefix] >= 1:
				del Orthogroup_cont[i][speclist_list.index(prefix)]
				if all(x == 0 for x in Orthogroup_cont[i]):
					unique_paralogs += Orthogroup_dict[i][prefix]
					do_Orthology.append(i)

	Other_orthology = list(set(Orthogroup_list).difference(set(do_Orthology)))
	
	for i in Other_orthology:
		Other_orthology_gene += Orthogroup_dict[i][prefix]


	unclass_gene = len(total_genes_list) - (single_copy_gene + multiple_copy_gene + unique_paralogs + Other_orthology_gene)



	print(f"single_copy_gene:{single_copy_gene}")
	print(f"multiple_copy_gene:{multiple_copy_gene}")
	print(f"unique_paralogs:{unique_paralogs}")
	print(f"Other_orthology_gene:{Other_orthology_gene}")
	print(f"unclass_gene:{unclass_gene}")
	#print((single_copy_gene + multiple_copy_gene + unique_paralogs + Other_orthology_gene))
	# print(len(total_genes_list))

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s GeneCount pep_fa prefix_in_GeneCount"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3])




