#获取对应基因家族的基因id

import re
import os
import sys

def main(Orthogroups_txt,fam_list):
	Orthogroups_dict = dict()
	with open(Orthogroups_txt) as inf:
		for lines in inf.readlines():
			line = lines.strip().split()
			# print(line[0])
			Orthogroups_dict[line[0].strip(":")]=line[1:]

	with open(fam_list) as inf:
		for lines in inf.readlines():
			line = lines.strip()
			print("\t".join(Orthogroups_dict[line]))

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s Orthogroups_txt fam_list"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])
