"""
1. 执行iqtree的sh文件, 为supertree 使用
"""

import re
import sys
import os 

iqtree = "path/to/home/anaconda3/envs/biosofeware/bin/iqtree"

def main(indir):
	indir = os.path.abspath(indir)
	fullfile_list = list()
	for file in os.listdir(indir):
		fullfile = os.path.join(indir,file)
		if fullfile.endswith("aln-cln"):
			fullfile_list.append(fullfile)

	for i in fullfile_list:
		print("{iqtree} -s {i} -m MFP -B 1000 -bnni -T AUTO -redo".format(iqtree=iqtree,i=i))

if __name__=="__main__":
	if len(sys.argv) != 2:
		print("!!! 修改@")
		print("python %s gblock_trimal_dir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])
