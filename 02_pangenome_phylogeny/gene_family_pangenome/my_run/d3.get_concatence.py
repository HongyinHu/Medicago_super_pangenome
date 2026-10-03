# 合并cds序列

import re
import os 
import sys
from Bio import SeqIO

def main(indir):
	seq_dict = dict()
	for dirpath, dirnames, files in os.walk(indir):
		for file in files:
			if file.endswith("aln.cds-gb"):
				fullfile = os.path.join(dirpath,file)

				for seq in SeqIO.parse(fullfile,"fasta"):
					seqid = seq.id.strip().split("|")[0]
					seq_seq = str(seq.seq).strip().replace(" ","")
					seq_dict.setdefault(seqid,[]).append(seq_seq)

	for k,v in seq_dict.items():
		res_seq = "".join(v)
		print(f">{k}\n{res_seq}")

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s indir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])