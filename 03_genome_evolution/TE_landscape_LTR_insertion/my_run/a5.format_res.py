#将结果进行特定格式格式化

import re
import os
import sys

species = "Msa_T2T","Mcr_468","Mma_457","Mtr_500","Mpr_410","Mpo_200","Mor_22","Mse_461","Mca_474","Mlu_395","Msu_472","Med_482","Mfi_46","Mra_436","Mru_300","Mar_100","Mla_454","Msa_ZM4"

def main(file):
	res_dict = dict()
	with open(file) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			res_dict[(line[0],line[1])] = line[2]

	for i in species:
		for ii in species:
			if (i,ii) in res_dict:
				print("\t".join([i,ii,res_dict[(i,ii)]]))
			elif (ii,i) in res_dict:
				print("\t".join([i,ii,res_dict[(ii,i)]]))
			else:
				print("\t".join([i,ii,"NA"]))

main(sys.argv[1])

