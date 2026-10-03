# 获取基因的功能结构域ID

import re
import os
import sys

def main(file):
	gene_domin = dict()

	with open(file) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			try:
				domin_ID = re.search(r"(IPR\d+)",lines).group(1)
				geneid = re.search(r"(Chr\d{5})",line[0]).group(1)
				gene_domin.setdefault(geneid,[]).append(domin_ID)
			except:
				pass

	for k,v in gene_domin.items():
		# v.insert(0,k)
		print(k+"\t"+",".join(v))

main(sys.argv[1])