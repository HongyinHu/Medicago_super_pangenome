# 将Orthogroups.GeneCount.tsv文件修改为cafe输入文件格式

import re
import os 
import sys

def main(GeneCount):
	with open(GeneCount) as inf:
		for index,lines in enumerate(inf.readlines()):
			line = lines.strip().split()
			line = line[0:-1]
			if index == "0":
				print("Descript"+"\t"+"\t".join(line))
			else:
				print("Null"+"\t"+"\t".join(line))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s GeneCount"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])



