# Orthogroup_Sequences中获取低拷贝同源基因集

import re
import os
import sys
from Bio import SeqIO

def main(data_dir,Ortho_seq_dir,outdir):
	species_list = list()
	for file in os.listdir(data_dir):
		species_name = file.strip(".fa")
		species_list.append(species_name)

	for file in os.listdir(Ortho_seq_dir):
		fullfile = os.path.join(Ortho_seq_dir,file)
		
		species_file = list()
		for seq in SeqIO.parse(fullfile,"fasta"):
			species_name = (seq.id).strip().split("|")[0]
			species_file.append(species_name)

		all_in_range = all( 1 <= species_file.count(i) <= 5 for i in species_list)
		if all_in_range:
			print(fullfile)
			os.system(f"cp {fullfile} {outdir}")

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s data_dir Ortho_seq_dir outdir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3])


