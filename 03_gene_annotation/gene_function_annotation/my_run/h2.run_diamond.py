# 使用diamond进行序列比对

import re
import os
import sys

script = "path/to/home/anaconda3/envs/biosofeware/bin/diamond"

def main(database,query,type_in="pep"):
	# query = os.path.abspath(query)

	if type_in == "pep":
		cmd = f"conda run -n biosofeware diamond blastp --threads 30 --db {database} --query {query} --out {query}_blastp.tab --outfmt 6 --sensitive --evalue 1e-5 --quiet"
		print(cmd)
	else:
		cmd = f"conda run -n biosofeware diamond blastx --threads 30 --db {database} --query {query} --out {query}_blastx.tab --outfmt 6 --sensitive --evalue 1e-5 --quiet"
		print(cmd)

def multip_run(indir,database,type_in):
	for dirpath,dirnames,files in os.walk(indir):
		for file in files:
			if file.endswith("pep.format"):
				fullfile = os.path.join(dirpath,file)
				main(database,fullfile,type_in)

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s indir database_no_prefile seqType[pep|cds]"%sys.argv[0])
		sys.exit(0)
	multip_run(sys.argv[1],sys.argv[2],sys.argv[3])