#对TE注释文件中的类别进行分类，获得相应的bed文件

import re
import os
import sys

def main(TE_gff):
	DNA_TE = ["DNA/DTA","DNA/DTC","DNA/DTH","DNA/DTM","DNA/DTT","MITE/DTA","MITE/DTC","MITE/DTH","MITE/DTM","MITE/DTT"]

	class_dict = dict()
	with open(TE_gff) as inf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split()
				type_TE = re.search(r"Classification=(.*?);",lines).group(1)
				if type_TE in DNA_TE:
					type_TE = "DNA_TEs"
				class_dict.setdefault(type_TE,[]).append((line[0],line[3],line[4],line[6]))

	for k,v in class_dict.items():
		filename = re.sub("/","_",k)
		with open(f"{filename}.TE.bed","w") as ouf:
			for i in v:
				print("\t".join(i),file=ouf)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s TE_gff"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])
