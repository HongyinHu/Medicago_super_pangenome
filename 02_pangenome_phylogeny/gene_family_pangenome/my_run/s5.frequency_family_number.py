# 计算基因家族包含物种数量频率分布

import re
import os 
import sys

def main(OrthoTsv):
	pan_geneFamily = dict()
	frequency_family = dict()

	with open(OrthoTsv) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index > 0:
				line = re.split(r"[\s,\t]+",lines)

				for i in line[1::]:
					if i:
						pan_geneFamily.setdefault(line[0].strip(),set()).add(i.strip().split("|")[0])

	for k,v in pan_geneFamily.items():
		frequency = len(v)
		frequency_family.setdefault(frequency,[]).append(k)

	for k,v in frequency_family.items():
		print(k,len(v))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s OrthoTsv_file"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])
	

