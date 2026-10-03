# 多序列比对并获取密码子序列

import re
import os 
import sys

script1 = "path/to/home/anaconda3/envs/biosofeware/bin/mafft"
script2 = "path/to/home/sofeware/PAL2NAL/pal2nal.pl"
script3 = "path/to/home/anaconda3/envs/biosofeware/bin/Gblocks"

def main(indir):

	for dirpath,dirnames,files in os.walk(indir):
		for file in files:
			if file.endswith("pep"):
				fullfile = os.path.join(dirpath,file)
				fullfile_cds = re.sub("fa.pep","fa.cds",fullfile)

				cmd1 = f"{script1} --auto {fullfile} > {fullfile}.aln.pep;"
				cmd2 = f"{script2} {fullfile}.aln.pep {fullfile_cds} -output fasta > {fullfile_cds}.aln.cds;"
				cmd3 = f"{script3} {fullfile_cds}.aln.cds -t c"
				print(" ".join([cmd1,cmd2,cmd3]))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s indir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])