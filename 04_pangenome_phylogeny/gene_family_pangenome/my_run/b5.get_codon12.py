#提取密码子的第1，第2位

import re
import os
import sys
from Bio import SeqIO

def main(indir):
	for file in os.listdir(indir):
		if file.endswith("pal2nal.fa"):
			fullfile = os.path.join(indir,file)
			with open(f"{fullfile}.codon12.fa","w") as ouf:
				for seq in SeqIO.parse(fullfile,"fasta"):
					sequence = str(seq.seq)
					res_seq = ""
					codon12 = dict()

					if len(sequence) % 3 != 0:
						print("序列长度不是3的倍数,无法提起")
					else:
						codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
						first_two_bases = [codon[:2] for codon in codons]

						for base in first_two_bases:
							res_seq += base

					print(f">{seq.id}\n{res_seq}",file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s indir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])