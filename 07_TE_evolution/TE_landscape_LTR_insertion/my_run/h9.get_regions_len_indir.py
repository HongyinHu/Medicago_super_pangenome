#获取区间长度累加之和
import re
import os
import sys
from Bio import SeqIO

def main(TE_indir,shared_bed):
	shared_len = 0
	with open(shared_bed) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			lens = int(line[2])-int(line[1])
			shared_len += lens

	all_TE_len = 0
	TE_dict = dict()
	for bed_file in os.listdir(TE_indir):
		if bed_file.endswith("merge.bed.shared.bed"):
			TE_len = 0
			with open(bed_file) as inf:
				for lines in inf.readlines():
					line = lines.strip().split()
					lens = int(line[2])-int(line[1])
					TE_len += lens

			all_TE_len += TE_len
			name = re.split(r"\.",bed_file)[0]
			TE_dict[name]=TE_len

	for k,v in TE_dict.items():
		print(k,v,"%.2f%%"%((v/all_TE_len)*100))



if __name__ == '__main__':
	if len(sys.argv) !=  3:
		print("python %s bed_file_indir shared_bed"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])
