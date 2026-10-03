# 将文件夹下的pep转为对应的cds序列

import re
import os
import sys
from Bio import SeqIO

def main(indir,cds_fa):
	cds_dict = dict()
	for seq in SeqIO.parse(cds_fa,"fasta"):
		cds_dict[seq.id] = seq.seq

	for file in os.listdir(indir):
		if file.endswith("align"):
			fullfile = os.path.join(indir,file)
			with open(f"{fullfile}.cds.fa","w") as ouf:
				for seq in SeqIO.parse(fullfile,"fasta"):
					print(f">{seq.id}\n{cds_dict[seq.id]}",file=ouf)

			cmd = f"pal2nal.pl {fullfile} {fullfile}.cds.fa -output fasta > {fullfile}.pal2nal.fa"
			# os.system(cmd)
			print(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s indir all_cds_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])

