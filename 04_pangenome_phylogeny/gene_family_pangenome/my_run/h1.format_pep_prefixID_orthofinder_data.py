# 为系统发育树候选物种蛋白序列添加物种名前缀

import re
import os 
import sys
from Bio import SeqIO


def main(indir,outdir):
	for file in os.listdir(indir):
		if file.endswith("fa"):
			prefix = re.split(r"\.",file)[0]
			print(prefix)
			fullfile = os.path.join(indir,file)

			with open("%s/%s.fa"%(outdir,prefix), "w") as ouf:
				for seq in SeqIO.parse(fullfile,"fasta"):
					if re.search(r"\.",str(seq.seq)):
						seq_seq = re.sub(r"\.","",str(seq.seq))
					else:
						seq_seq = seq.seq

					if re.search(r"\|", seq.id):
						newid = seq.id 
					else:
						newid = "|".join([prefix, seq.id])
					
					print(">"+newid+"\n"+seq_seq,file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s indir outdir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2])
