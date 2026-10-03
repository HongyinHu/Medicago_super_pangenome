# 获取蛋白序列的blastp比对结果

import re
import os 
import sys

script1 = "path/to/home/sofeware/CEGMA_v2/depend_sofe/blast+/ncbi-blast-2.10.0+/bin/makeblastdb"
script2 = "path/to/home/sofeware/CEGMA_v2/depend_sofe/blast+/ncbi-blast-2.10.0+/bin/blastp"

def main(pep_A, pep_B):
	pep_A = os.path.abspath(pep_A)
	pep_B = os.path.abspath(pep_B)

	cmd1 = f"diamond makedb --in {pep_B} -d {os.path.basename(pep_B)}"
	cmd2 = f"diamond blastp --threads 10 --db {os.path.basename(pep_B)} --query {pep_A} --out all.blastp.txt --outfmt 6 --sensitive --max-target-seqs 10 --evalue 1e-5"
	print(cmd1)
	print(cmd2)


if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s pep_A pep_B_anc"%sys.argv[0])
		sys.exit(0)

	main(sys.argv[1], sys.argv[2])
