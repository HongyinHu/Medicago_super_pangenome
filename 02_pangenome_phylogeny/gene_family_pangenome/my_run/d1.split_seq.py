# 对每个单拷贝基因家族进行处理，创建该家族单独文件夹，包含pep和cds序列

import re
import os 
import sys
import shutil
from Bio import SeqIO

def main(singlecopy_indir,all_cds):
	
	for file in os.listdir(singlecopy_indir):
		dirname = file.strip(".fa")
		if not os.path.exists(f"align/{dirname}"):
			os.makedirs(f"align/{dirname}")

		seq_id_list = list()
		fullfile = os.path.join(singlecopy_indir,file)
		shutil.copy(fullfile,f"align/{dirname}/{file}.pep")
		for seq in SeqIO.parse(fullfile,"fasta"):
			seq_id_list.append(seq.id)

		with open(f"align/{dirname}/{file}.cds","w") as ouf:
			seq_dict = dict()
			for seq in SeqIO.parse(all_cds,"fasta"):
				seq_dict[seq.id] = seq.seq

			for i in seq_id_list:
				print(f">{i}\n{seq_dict[i]}",file=ouf)
if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s singlecopy_indir all_cds"%sys.argv[0])
		sys.exit(0)

	main(sys.argv[1],sys.argv[2])