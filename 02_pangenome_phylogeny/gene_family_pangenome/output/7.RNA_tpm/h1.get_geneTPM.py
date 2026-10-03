"""
1. 从stringtie ballgown 中产生的样品注释文件中提取转录本的表达量
"""

import re
import os 
import sys

def main(indir):
	indir = os.path.abspath(indir)

	file_list = list()
	file_name_list = list()
	gene_list = list()
	sample_TPM = dict()

	for file in os.listdir(indir):
		if re.search(r"gtf",file):
			file_name = re.search(r"(.*?).RNA.gtf",file).group(0)
			full_file = os.path.join(indir,file)
			file_name_list.append(file_name)
			file_list.append(file)
			
			with open(full_file) as inf, open(f"{full_file}.tpm","w") as ouf:
				for lines in inf.readlines():
					if lines.startswith("#"): continue
					line = lines.strip().split("\t")
					if line[2] == "transcript":
						transcript_id = re.search(r"transcript_id\s+\"(.*?)\";",line[-1]).group(1)
						TPM = re.search(r"TPM\s+\"(.*?)\";",line[-1]).group(1)
						print("\t".join([transcript_id,TPM]),file=ouf)
						

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("python %s gtf_dir"%sys.argv[0])
		sys.exit(0)

	main(sys.argv[1])