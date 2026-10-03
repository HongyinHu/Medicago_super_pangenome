# 产生两两组合的基因对，用于计算成对的Ka/KS
# s9 废掉就可
# perl path/to/home/sofeware/ParaAT2.0/ParaAT.pl -h ./res_homologs.txt -n ./res_cds.fa -a ./res_pep.fa  -p proc -m mafft -f axt -g -k -o result_dir

import os
import sys
import itertools
import random
import re

script = "perl path/to/home/sofeware/ParaAT2.0/ParaAT.pl"

def gene_combinations(gene_list,num_pairs):
	all_pairs = list(itertools.combinations(gene_list,2))
	if len(all_pairs) >= int(num_pairs):
		sample_pairs = random.sample(all_pairs,int(num_pairs))
		return sample_pairs


def generate_run_cmd(homologous_list_file,cds_fa,pep_fa,random_num):
	if not os.path.exists("homo_split"):
		os.mkdir("homo_split")

	# 两两排列组合
	num = 1
	homologous_list = list()
	with open(homologous_list_file) as inf:
		for lines in inf.readlines():
			line = lines.strip().split("\t")
			random_pairs = gene_combinations(line,random_num)
			if random_pairs:
				with open(f"homo_split/genefam_{num}.txt","w") as ouf:
					for i in random_pairs:
						print("\t".join(i),file=ouf)

			num += 1

	for file in os.listdir("homo_split"):
		if file.endswith("txt"):
			fullfile = os.path.join("homo_split",file)
			out_dir = "homo_res/%s"%re.sub(".txt","",file)
			# out_dir = "homo_res"
			
			if not os.path.exists("homo_res"):
				os.mkdir("homo_res")

			cmd = f"{script} -h {fullfile} -n {cds_fa} -a {pep_fa} -m mafft -p proc -f axt -k -o {out_dir}"
			print(cmd)


if __name__ == '__main__':
	if len(sys.argv) != 5:
		print("python %s homologous_list_file cds_fa pep_fa random_sample_num"%sys.argv[0])
		sys.exit(0)
	generate_run_cmd(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])