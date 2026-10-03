#获取在共享区域的LTR的序列，并集有一段在即可

import re
import os
import sys
from Bio import SeqIO

def main(genome_fa,lLTR_bed,rLTR_bed,prefix):
	genome_dict = dict()
	for seq in SeqIO.parse(genome_fa,"fasta"):
		genome_dict[seq.id]=str(seq.seq)

	

	all_LTR = dict()

	with open(lLTR_bed) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			name_l = line[-2].lstrip("l")
			all_LTR.setdefault(name_l,{})[line[-2]]=(line[0],line[1],line[2])

	with open(rLTR_bed) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			name_l = line[-2].lstrip("r")
			all_LTR.setdefault(name_l,{})[line[-2]]=(line[0],line[1],line[2])

	if not os.path.exists(f"LTR_{prefix}"):
		os.mkdir(f"LTR_{prefix}")
	
	for i in all_LTR:
		with open(f"LTR_{prefix}/{i}.fa","w") as ouf:
			ltr_l_coords = all_LTR[i][f"l{i}"]
			ltr_r_coords = all_LTR[i][f"r{i}"]

			ltr_l_seq = genome_dict[ltr_l_coords[0]][int(ltr_l_coords[1]):int(ltr_l_coords[2])]
			ltr_r_seq = genome_dict[ltr_r_coords[0]][int(ltr_r_coords[1]):int(ltr_r_coords[2])]

			print(f">l{i}\n{ltr_l_seq}",file=ouf)
			print(f">r{i}\n{ltr_r_seq}",file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 5:
		print("python %s genome_fa lLTR_bed rLTR_bed prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

			




