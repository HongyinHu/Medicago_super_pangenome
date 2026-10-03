# 获取各个数据库中注释到基因数目

import re
import os
import sys

def main(indir,name_database):
	for dirpath, dirnames, files in os.walk(indir):
		for file in files:
			if file.endswith("format.interproscan.out"):
				fullfile = os.path.join(dirpath,file)
				name = file.strip().split(".")[0]
				cmd = f"cut -f 1 {fullfile} | sort | uniq > {dirpath}/{name}.{name_database}.genelist"
				print(cmd)
				os.system(cmd)

def main2(indir,name_database):
	for dirpath, dirnames, files in os.walk(indir):
		for file in files:
			if file.endswith("format.pfam.out"):
				fullfile = os.path.join(dirpath,file)
				name = file.strip().split(".")[0]
				with open(fullfile) as inf, open(f"{dirpath}/{name}.{name_database}.genelist","w") as ouf:
					gene_set = set()
					for lines in inf.readlines():
						if not lines.startswith("#"):
							line = lines.strip().split()
							if len(line) > 4:
								gene_set.add(line[0])
					for i in gene_set:
						print(i,file=ouf)

				# cmd = f"cut -f 1 {fullfile} | sort | uniq > {dirpath}/{name}.{name_database}.genelist"
				# print(cmd)
				# os.system(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s indir name_database"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])