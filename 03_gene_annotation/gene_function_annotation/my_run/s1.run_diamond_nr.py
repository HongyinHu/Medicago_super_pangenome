# 使用diamond 对蛋白功能进行功能注释

import re
import os 
import sys

script = "path/to/home/anaconda3/envs/biosofeware/bin/diamond"

def main(diamond_dmnd, pep_fa):
	cmd = f"{script} blastp -d {diamond_dmnd} -q {pep_fa} -o res.blastp -p 20 -f 6 -e 0.00001"
	print(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s diamond_dmnd pep_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2])