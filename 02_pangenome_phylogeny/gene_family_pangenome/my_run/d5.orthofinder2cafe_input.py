# 将Orthogroups.GeneCount.tsv整理成cafe的输入格式

import re
import os 
import sys

def main(file):
	with open(file) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index == 0:
				header = lines.strip().split()[0:-1]
				header.insert(0,"Descript")
				print("\t".join(header))
			if index > 0:
				line = lines.strip().split()
				need_line = line[0:-1]
				need_line.insert(0,"Null")
				print("\t".join(need_line))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s Orthogroups.GeneCount.tsv"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])
