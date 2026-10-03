#获取beb文件交集，如果A文件完全被B文件覆盖，则输出A文件坐标

import re
import os
import sys

def main(a_bed,b_bed):
	B_bed_list = list()
	with open(b_bed) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			B_bed_list.append(line)

	with open(a_bed) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			for i in B_bed_list:
				if line[0] == i[0] and line[1] >= i[1] and line[2] <= i[2]:
					print(lines,end="")

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s a_bed b_bed"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])