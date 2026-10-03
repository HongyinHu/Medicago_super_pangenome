# 获取shared TE 区域
import re
import os
import sys

def main(TE_bed_dir, homology_bed):
	homology_dict = dict()
	with open(homology_bed) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			homology_dict.setdefault(line[0],[]).append([line[1],line[2]])

	for TE_bed in os.listdir(TE_bed_dir):
		if TE_bed.endswith("TE.convert.bed"):

			all_TE = list()
			share_TE = list()
			with open(TE_bed) as inf:
				for lines in inf.readlines():
					all_TE.append(lines)
					
					line = lines.strip().split()
					for i in homology_dict[line[0]]:
							if line[1] >= i[0] and line[2] <= i[1]:
								share_TE.append(lines)

			specific_list = set(all_TE) - set(share_TE)
			
			with open(f"{TE_bed}.shared.bed","w") as ouf1, open(f"{TE_bed}.specific.bed","w") as ouf2:
				for i in specific_list:
					print(i,end="",file=ouf2)

				for ii in share_TE:
					print(ii,end="",file=ouf1)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s TE_bed_dir homology_bed"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])



			


