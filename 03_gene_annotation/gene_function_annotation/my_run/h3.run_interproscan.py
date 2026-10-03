# 使用interproscan进行GO注释

import re
import os
import sys

#script = "path/to/software/interproscan/interproscan-5.30-69.0/interproscan.sh"
script = "path/to/home/sofeware/my_interproscan/interproscan-5.30-69.0/interproscan.sh"

def main(query_pep,outdir):
	cmd = f"conda run -n java {script} -f tsv -i {query_pep} -cpu 20 --highmem -o {query_pep}.interproscan.out -iprlookup -goterms -pa -td {outdir}/temp"
	print(cmd)

def multip_run(indir):
	for dirpath,dirnames,files in os.walk(indir):
		for file in files:
			if file.endswith("pep.format"):
				fullfile = os.path.join(dirpath,file)
				main(fullfile,os.path.dirname(fullfile))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s indir"%sys.argv[0])
		sys.exit(0)
	multip_run(sys.argv[1])
