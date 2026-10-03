# 运行busco评估基因组组装和蛋白注释

import re
import os 
import sys

database = "path/to/home/busco_version/embryophyta_1614_odb10"
# busco = "path/to/home/anaconda3/envs/busco/bin/busco"
busco = "path/to/home/anaconda3/envs/biosofeware/bin/busco"

def main(seq_type, fasta):
	fasta = os.path.abspath(fasta)
	# print("conda activate busco")
	print("{busco} -i {fasta} -l {database} -o busco_assess --out_path ./ -c 30 -m {seq_type} -q --offline -f".format(busco=busco, fasta=fasta, 
		database=database, seq_type=seq_type))

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s mode[genome|proteins] fasta"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2])

