#bed文件进行sort和merge操作

import re
import os
import sys

def main(bed_indir):
	for bed_file in os.listdir(bed_indir):
		print(bed_file)
		if bed_file.endswith("TE.bed"):
			cmd0 = f"cut -f 1,2,3 {bed_file} > {bed_file}.cut"
			cmd1 = f"sort -k1,1 -k2,2n {bed_file}.cut > {bed_file}.sort.bed"
			cmd2 = f"bedtools merge -i {bed_file}.sort.bed > {bed_file}.sort.merge.bed"
			os.system(cmd0)
			os.system(cmd1)
			os.system(cmd2)
			print("h7 run succeed!")

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s bed_indir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])