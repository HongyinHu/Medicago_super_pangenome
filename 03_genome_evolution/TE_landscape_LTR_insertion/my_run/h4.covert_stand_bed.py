#转换链信息，将负链的位置坐标转到相应的正链位置坐标上
import re
import os
import sys
from Bio import SeqIO

def main(bed_file_dir, genome_fa):
	chr_len = dict()
	for seq in SeqIO.parse(genome_fa,"fasta"):
		chr_len[seq.id] = len(seq.seq)

	for file in os.listdir(bed_file_dir):
		if file.endswith("bed"):
			newfile = re.sub("bed","convert.bed",file)
			with open(file) as inf, open(newfile,"w") as ouf:
				for lines in inf:
					line = lines.strip().split()
					if line[-1] == "-":
						new_start = chr_len[line[0]] - int(line[2])
						new_end = chr_len[line[0]] - int(line[1])
						new_strand = "+"
						print("\t".join([line[0],str(new_start),str(new_end),new_strand]),file=ouf)
					else:
						print("\t".join(line),file=ouf)


if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s bed_file_dir genome_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2])


