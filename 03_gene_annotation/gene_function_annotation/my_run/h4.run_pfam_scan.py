# 使用pfam_scan搜索pfam数据库

import re
import os
import sys

script = "path/to/home/anaconda3/envs/pfam/bin/pfam_scan.pl"

def main(query_pep,pfam_dir):
	cmd = f"conda run -n pfam {script} -fasta {query_pep} -dir {pfam_dir} -outfile {query_pep}.pfam.out -cpu 30"
	print(cmd)

def multip_run(indir,database):
	for dirpath,dirnames,files in os.walk(indir):
		for file in files:
			if file.endswith("pep.format"):
				fullfile = os.path.join(dirpath,file)
				main(fullfile, database)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s indir pfam_dir"%sys.argv[0])
		sys.exit(0)
	multip_run(sys.argv[1],sys.argv[2])