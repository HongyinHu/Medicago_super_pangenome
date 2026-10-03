# 对每个基因组蛋白序列运行interproscan进行interpro结构域注释

import re
import os 
import sys
from Bio import SeqIO
import time

def split_fa(pep_fa):
	file_num = 1
	seq_list = list()

	if not os.path.exists("split_fa"):
		os.mkdir("split_fa")

	for seq in SeqIO.parse(pep_fa,"fasta"):
		if len(seq_list) <= 100:
			seq_list.append(seq)

		else:
			SeqIO.write(seq_list,"split_fa/split_pep_%s.fa"%file_num,"fasta")
			file_num += 1
			seq_list = []

def main(pep_fa):
	split_fa(pep_fa)

	if not os.path.exists("interproscan_out"):
		os.mkdir("interproscan_out")

	with open("run_interproscan.sh","w") as ouf:
		for file in os.listdir("split_fa"):
			fullfile = os.path.join("split_fa",file)

			cmd = "interproscan.sh -appl TIGRFAM,ProDom,Hamap,SMART,ProSiteProfiles,ProSitePatterns,SUPERFAMILY,PRINTS,Gene3D,PIRSF,Pfam,Coils \
-f tsv -t p -i %s  -o interproscan_out/interpro_anno_%s.out  -iprlookup -goterms -pa -td temp && echo '** go annotation done %s **' >> run_interproscan.log"%(fullfile,file,file)
			
			print(cmd,file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s pep_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1]) 
