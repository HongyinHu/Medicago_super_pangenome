# 使用merqury评估基因组组织质量

import re
import os 
import sys
from Bio import SeqIO

script1 = "path/to/home/sofeware/merqury/best_k.sh"
script2 = "path/to/home/anaconda3/envs/biosofeware/bin/meryl"
script3 = "path/to/home/sofeware/merqury/merqury.sh"

def main(genome_fa,paired_read1,paired_read2):
	# 确定合适的kmer大小
	genome_size = 0
	for seq in SeqIO.parse(genome_fa,"fasta"):
		genome_size += len(str(seq.seq))

	cmd1 = f"{script1} {genome_size}"
	k = round(float(os.popen(cmd1).read().split("\n")[2]))
	
	#构建kmer的db
	cmd2 = f"{script2} k={k} count output {paired_read1}.meryl {paired_read1}"
	cmd3 = f"{script2} k={k} count output {paired_read2}.meryl {paired_read2}"
	cmd4 = f"{script2} union-sum output read-db.meryl {paired_read1}.meryl {paired_read2}.meryl"
	cmd5 = f"{script3} read-db.meryl {genome_fa} out"
	print("\n".join([cmd2,cmd3,cmd4,cmd5]))

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s assembly_fa paired_read1 paired_read2"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2], sys.argv[3])

