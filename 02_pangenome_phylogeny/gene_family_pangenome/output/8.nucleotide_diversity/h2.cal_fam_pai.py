# 计算泛基因组家族类型的核酸多样性

import re
import os
import sys

cal_script = "path/to/project/1.orthology_family_2/output/8.nucleotide_diversity/nucleotide_diversity.R"

def main(msa_dir):
	file_para = list()
	for file in os.listdir(msa_dir):
		fullfile = os.path.join(msa_dir,file)
		with open(fullfile) as inf:
			for index,lines in enumerate(inf.readlines()):
				if index == 0:
					sample_num = lines.strip().split()[0]
					cmd=f"Rscript {cal_script} {fullfile} {sample_num} {fullfile}.nuc.div.txt"
					print(cmd)
				
if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s msa_indir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])