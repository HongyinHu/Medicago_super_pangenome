# 修改终止符号

import re
import os 
import sys
from Bio import SeqIO

def main(pep_fa):
	for seq in SeqIO.parse(pep_fa,"fasta"):
		newseq = re.sub(r"\.","",str(seq.seq))
		print(f">{seq.id}\n{newseq}")

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s pep_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])