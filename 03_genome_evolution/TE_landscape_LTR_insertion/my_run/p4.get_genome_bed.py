#获取基因组bed文件，为差异就是非共享区域

import re
import os
import sys
from Bio import SeqIO

def main(fasta):
	for seq in SeqIO.parse(fasta,"fasta"):
		print("\t".join([seq.id,"0",str(len(seq.seq)-1),"+"]))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s fasta"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])