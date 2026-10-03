# 获取系统发育树使用的单拷贝基因家族
#  -os Stop after writing sequence files for orthogroups

import re
import os 
import sys

script = "path/to/home/anaconda3/envs/biosofeware/bin/orthofinder"

def main(genome_peps_indir, out_dir):
	print(f"{script} -f {genome_peps_indir} -o {out_dir} -M msa -os -t 30")

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s genome_peps_indir out_dir"%sys.argv[0])
		sys.exit()
	main(sys.argv[1],sys.argv[2])