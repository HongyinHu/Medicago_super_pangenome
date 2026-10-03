# 过滤Mpr_Mpr_block_information.csv大小
import re
import os
import sys

def main(file):
	with open(file) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index > 0:
				line = lines.strip().split(",")
				if 0.5 <= float(line[10]) <= 1.5:
					print(lines,end="")
			else:
				print(lines,end="")

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s block_information_csv"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])
