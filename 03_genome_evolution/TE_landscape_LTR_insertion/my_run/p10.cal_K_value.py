#使用MEGAcc计算K值(sequence divergence rate between its 5` and 3` LTRs)

import re
import os
import sys

MEGA = "path/to/home/sofeware/megacc_11.0.13_i386/megacc"
align_mao = "path/to/project/17.TE_landscape/my_run/muscle_align_nucleotide.mao"
estimation_mao = "path/to/project/17.TE_landscape/my_run/distance_estimation_pairwise_nucleotide.mao"

def main(ltr_indir,step):

	if step == "1":
		if not os.path.exists(f"{ltr_indir}_res1"):
			os.mkdir(f"{ltr_indir}_res1")

		with open(f"{ltr_indir}.cmd1.sh","w") as ouf:
			for file in os.listdir(ltr_indir):
				if file.endswith(".fa"):
					fullfile = os.path.join(ltr_indir,file)
					cmd1 = f"{MEGA} -a {align_mao} -d {fullfile} -o {ltr_indir}_res1"
					# os.system(cmd1)
					print(cmd1,file=ouf)

	if step == "2":
		if not os.path.exists(f"{ltr_indir}_res2"):
			os.mkdir(f"{ltr_indir}_res2")

		with open(f"{ltr_indir}.cmd2.sh","w") as ouf1:
			for file in os.listdir(f"{ltr_indir}_res1"):
				if file.endswith(".meg"):
					fullfile = os.path.join(f"{ltr_indir}_res1",file)
					cmd2 = f"{MEGA} -a {estimation_mao} -d {fullfile} -o {ltr_indir}_res2"
					# os.system(cmd2)
					print(cmd2,file=ouf1)


if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s ltr_indir step"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2])




