#批量运行vmatch用于大规模聚类计算
# time: 20241124
# author: hhy

import re
import os
import sys
import shutil

script1 = "path/to/home/sofeware/vmatch-2.3.1-Linux_x86_64-64bit/mkvtree"
script2 = "path/to/home/sofeware/vmatch-2.3.1-Linux_x86_64-64bit/vmatch"

def main(indir):
	indir = os.path.abspath(indir)

	for dirpath,dianames,files in os.walk(indir):
		original_dir = os.getcwd()

		for file in files:
			if file == "ltr_juction.fa":
				os.chdir(dirpath)

				with open("run.sh","w") as ouf:
					if os.path.exists("70"): shutil.rmtree("70")
					os.mkdir(f"70")
					if os.path.exists("80"): shutil.rmtree("80")
					os.mkdir(f"80")
					if os.path.exists("90"): shutil.rmtree("90")
					os.mkdir(f"90")

					cmd1 = f"{script1} -db {file} -dna -pl -allout -v"
					cmd2 = f"{script2} -dbcluster 70 70 70/Cluster -s -identity 70 -exdrop 6 -seedlength 10 -d {file}"
					cmd3 = f"{script2} -dbcluster 80 80 80/Cluster -s -identity 80 -exdrop 5 -seedlength 15 -d {file}"
					cmd4 = f"{script2} -dbcluster 90 90 90/Cluster -s -identity 90 -exdrop 4 -seedlength 20 -d {file}"

					print(" && ".join([cmd1,cmd2,cmd3,cmd4]),file=ouf)
					# print(" && ".join([cmd1,cmd4]),file=ouf)

				os.chdir(original_dir)

main(sys.argv[1])