# 计算在随意组合情况下，pan和core的基因家族数量变化趋势

import re
import os 
import sys
import itertools
from multiprocessing import Pool

def cal_combind(combind):
	core_list = list()
	pan_list = list()
	
	for k,v in pan_geneFamily.items():
		if all(item in v for item in combind):
			core_list.append(k)
		
		if any(item in v for item in combind):
			pan_list.append(k)
	
	print(len(combind),len(pan_list),len(core_list))

def main(OrthoTsv,Core_file):
	global pan_geneFamily
	pan_geneFamily = dict()
	
	species_name = set()
	
	global Core_list
	Core_list = list()

	with open(OrthoTsv) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index > 0:
				line = re.split(r"[\s,\t]+",lines)

				for i in line[1::]:
					if i:
						pan_geneFamily.setdefault(line[0].strip(),set()).add(i.strip().split("|")[0])
						species_name.add(i.strip().split("|")[0])

	with open(Core_file) as inf:
		for lines in inf.readlines():
			Core_list.append(lines.strip())

	pool = Pool(50)
	# for i in range(2, 19):
	i = 18
	combination = list(itertools.combinations(species_name,i))
	pool.map(cal_combind,combination)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s Orthogroup_Tsv pangenome_Core"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])

