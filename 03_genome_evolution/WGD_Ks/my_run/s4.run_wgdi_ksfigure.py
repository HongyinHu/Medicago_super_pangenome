# 书写wgdi配置文件和运行各个命令

import re
import os 
import sys

def run_ksfigure(all_ks,prefix):
	conf9 = """
[ksfigure]
ksfit = {all_ks}
labelfontsize = 15
legendfontsize = 15
xlabel = none
ylabel = none
title = none
area = 0,5
figsize = 10,10
shadow = false
savefig = {prefix}.ksfigure.pdf
""".format(all_ks=all_ks, prefix=prefix)
	cmd9 = "wgdi -kf total.conf"
	return (conf9, cmd9)

def main(all_ks, genome1_name, genome2_name):
	all_ks = os.path.abspath(all_ks)
	prefix = "_".join([genome1_name, genome2_name])
	a9 = run_ksfigure(all_ks,prefix)
	with open("total_ksfigure.conf","w") as ouf:
		print(a9[0],file=ouf)

	print(a9[1])

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("python %s all_ks genome1_name genome2_name"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2], sys.argv[3])


