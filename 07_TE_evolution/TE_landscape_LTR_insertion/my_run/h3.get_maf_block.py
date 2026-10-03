# 获取多个基因组序列比对blocks, 每个block至少共享在两个物种中，最后使用bedtools进行合并处理

import re
import os
import sys

def search_block(i,species_name):
	res_list = list()
	if re.search(species_name,i) and len(i.split("\n")) >= 3:
		for lines in i.split("\n"):
			if lines.startswith("s") and re.search(species_name, lines):
				line = lines.strip().split()
				# print(line)
				Chr = re.split(r"\.",line[1])[1]
				END = int(line[2]) + int(line[3])
				strand = line[4]
				res_list.append("\t".join([Chr,line[2],str(END),strand]))

	return res_list

def main(maf_file):
	names = ["Mar_100","Mru_300","Mla_454","Mca_474","Mcr_468","Msa_T2T","Msa_zm4","Mma_457","Mpr_410","Mtr_500","Mpo_200","Mor_22","Mse_461","Mlu_395","Msu_472","Med_482","Mfi_46","Mra_436"]
	with open(maf_file) as inf:
		blocks = inf.read().split("\n\n")
		for species_name in names:
			with open(f"{species_name}.blocks.bed","w") as ouf:
				for i in blocks:
					res = search_block(i,species_name)
					for res_out in res:
						print(res_out,file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s maf_file"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])
