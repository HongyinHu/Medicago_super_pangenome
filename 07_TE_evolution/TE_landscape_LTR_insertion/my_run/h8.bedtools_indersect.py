#bed文件进行sort和merge操作

import re
import os
import sys

def main(bed_indir,shared_BED):
	for bed_file in os.listdir(bed_indir):
		if bed_file.endswith("sort.merge.bed"):
			cmd0 = f"bedtools intersect -a {shared_BED} -b {bed_file} > {bed_file}.shared.bed"
			os.system(cmd0)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s bed_indir shared_BED"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])