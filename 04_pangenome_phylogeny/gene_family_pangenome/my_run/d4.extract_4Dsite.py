# 提取四重简并为点从fasta文件中

import re 
import os
import sys
from Bio import SeqIO

def main(file):
	codon_table={
	'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
	'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
	'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
	'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
	'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
	'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
	'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
	'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
	
	wrt = []
	for rec in SeqIO.parse(file,"fasta"):
		i=0
		wrt.append('>'+str(rec.id)+'\n')
		
		while i < len(rec.seq):
			if str(rec.seq[i:i+3]) in codon_table:
				wrt.append(str(rec.seq[i+2]))
			i+=3
		wrt.append('\n')
	
	with open(file +'.4Dsite','w') as f:
		f.write(''.join(wrt))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s file"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])
