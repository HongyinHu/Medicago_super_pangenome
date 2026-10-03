# 获取pfam注视到的基因数量

import re
import os
import sys

def main(file):
	gene_set = set()

	with open(file) as inf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split()
				if len(line) > 3:
					gene_set.add(line[0])

	print(gene_set)
	print(len(gene_set))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s pfam_out"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])