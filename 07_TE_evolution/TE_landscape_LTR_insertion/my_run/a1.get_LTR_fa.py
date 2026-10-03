#获取fl_LTR序列
import re
import os
import sys
from Bio import SeqIO

def main(genome_fa,LTR_gff,prefix):
	genome_fa_dict = dict()
	for seq in SeqIO.parse(genome_fa,"fasta"):
		genome_fa_dict[seq.id] = str(seq.seq)

	LTR_dict = dict()
	with open(LTR_gff) as inf:
		for lines in inf.readlines():
			if re.search("long_terminal_repeat",lines):
				line = lines.strip().split()
				ids = re.search(r"ID=(.*?);",lines).group(1)
				ids_name = re.search(r"(LTR_\d+)",ids).group(1)
				if re.search("lLTR",ids):
					LTR_dict.setdefault(ids_name,{})["lLTR"] = [line[0], int(line[3])-1, int(line[3])+100,line[6]]
				if re.search("rLTR",ids):
					LTR_dict.setdefault(ids_name,{})["rLTR"] = [line[0], int(line[4])-100, int(line[4])+1,line[6]]

	# with open(f"{prefix}.lLTR.fa","w") as ouf1, open(f"{prefix}.rLTR.fa","w") as ouf2:
	with open(f"{prefix}.LTR.fa","w") as ouf:
		for k,v in LTR_dict.items():
			seql = genome_fa_dict[v["lLTR"][0]][v["lLTR"][1]:v["lLTR"][2]]
			seqr = genome_fa_dict[v["rLTR"][0]][v["rLTR"][1]:v["rLTR"][2]]
			print(f">{prefix}_{k}\n{seql}{seqr}",file=ouf)
			# print(f">{prefix}_r{k}\n{seqr}",file=ouf2)


if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s genome_fa LTR_gff prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2],sys.argv[3])


