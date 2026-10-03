# 获取单拷贝基因家族序列后，进行系统发育建树

import re
import os 
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq 

script1 = "path/to/home/anaconda3/envs/biosofeware/bin/mafft"
script2 = "path/to/home/anaconda3/envs/biosofeware/bin/Gblocks"
script3 = "path/to/home/anaconda3/envs/biosofeware/bin/trimal"
script4 = "path/to/home/sofeware/iqtree-1.6.12-Linux/bin/iqtree"

def main(single_copy_indir):
	if not os.path.exists("run_mafft"):
		os.mkdir("run_mafft")

	with open("r1.run_mafft.sh","w") as ouf1:
		for file in os.listdir(single_copy_indir):
			if file.endswith("fa"):
				fullfile = os.path.join(single_copy_indir,file)
				print("{script1} {fullfile} > run_mafft/{file}.align".format(script1=script1, fullfile=fullfile, file=file),file=ouf1)

	with open("r2.run_gblock.sh","w") as ouf2:
		for file in os.listdir("run_mafft"):
			if file.endswith("align"):
				fullfile = os.path.join("run_mafft", file)
				print("{script2} {fullfile} -b4=5 -b5=h".format(script2=script2, fullfile=fullfile),file=ouf2)

	with open("r3.run_trimal.sh","w") as ouf3:
		for file in os.listdir("run_mafft"):
			if file.endswith("align-gb"):
				fullfile = os.path.join("run_mafft",file)
				print("{script3} -in {fullfile} -out {fullfile}.aln-cln -gt 0.8".format(script3=script3, fullfile=fullfile),file=ouf3)

	all_species = set()
	species_dict = dict()
	
	for file in os.listdir("run_mafft"):
		if file.endswith("aln-cln"):
			fullfile = os.path.join("run_mafft",file)
			for seq in SeqIO.parse(fullfile, "fasta"):
				seqid = seq.id.strip().split("|")[0]
				all_species.add(seqid)

	for species in all_species:
		if species not in species_dict:
			species_dict[species] = []

	for file in os.listdir("run_mafft"):
		if file.endswith("aln-cln"):
			file_dict = dict()
			seq_len = 0

			fullfile = os.path.join("run_mafft",file)
			for seq in SeqIO.parse(fullfile, "fasta"):
				seqid = seq.id.strip().split("|")[0]
				file_dict[seqid] = str(seq.seq)
				seq_len = len(str(seq.seq))

			for k in all_species:
				if k in file_dict:
					species_dict[k].append(file_dict[k])
				else:
					seq_seq = "-"*seq_len
					species_dict[k].append(seq_seq)

	with open("concatence_iqtree.fa","w") as ouf4:
		for k,v in species_dict.items():
			print(">"+k+"\n"+"".join(v),file=ouf4)

	if os.path.exists("concatence_iqtree.fa"):
		with open("r4.run_iqtree.sh","w") as ouf5:
			print("{script4} -s concatence_iqtree.fa -m MFP -bb 1000 -bnni -redo -nt 50 -pre concatence_iqtree".format(script4=script4),file=ouf5)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s single_copy_indir_from_orthofinder"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])











