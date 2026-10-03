# 为WGDI输入文件过滤部分非必要染色体片段

import re
import os 
import sys

def main(file):
	with open(file) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			if re.search("random",line[0]):
				continue
			elif re.search("chrUn",line[0]):
				continue
			else:
				print(lines,end="")

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s file"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])
