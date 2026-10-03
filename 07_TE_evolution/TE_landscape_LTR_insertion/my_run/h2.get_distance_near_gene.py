#获取各类型TE和gene的上下游距离
import re
import os
import sys

def main(file):
	with open(file) as inf, open("all.TE.uptream.out","w") as ouf1, open("all.TE.downtream.out","w") as ouf2:
		for lines in inf.readlines():
			line = lines.strip().split()
			TE_type = re.search(r"([a-zA-Z_]+)\d+",line[3]).group(1)
			if 0 < abs(int(line[-1])) <= 10000:
				if int(line[-1]) > 0:
					print("\t".join([TE_type, str(int(line[-1])/1000)]),file=ouf2)
				else:
					print("\t".join([TE_type, str(int(line[-1])/1000)]),file=ouf1)
					
if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s all.TE.bed"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])