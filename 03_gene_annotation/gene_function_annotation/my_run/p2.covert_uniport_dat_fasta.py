# 将uniport dat格式转为fasta格式

import re
import os
import sys
from Bio import SeqIO

def main(dat_file):
	count = SeqIO.convert(dat_file,"swiss",dat_file+".fa","fasta")
	print("Coerted %i records"%count)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s swiss_dat"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])