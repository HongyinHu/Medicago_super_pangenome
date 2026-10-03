#获取各类型TE的比例

import re
import os
import sys

h6 = "path/to/project/17.TE_landscape/my_run/h6.covert_stand_bed_indir.py"
h7 = "path/to/project/17.TE_landscape/my_run/h7.bed_sort_merge.py"
h8 = "path/to/project/17.TE_landscape/my_run/h8.bedtools_indersect.py"
h9 = "path/to/project/17.TE_landscape/my_run/h9.get_regions_len_indir.py"

def main(prefix):
	cmd1 = f"python {h6} ./ ../{prefix}.genome.fa"
	cmd2 = f"python {h7} ./"
	cmd3 = f"python {h8} ./ ../{prefix}.shared.convert.sort.specific.bed"
	cmd4 = f"python {h9} ./ ../{prefix}.shared.convert.sort.specific.bed"

	os.system(cmd1)
	os.system(cmd2)
	os.system(cmd3)
	os.system(cmd4)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s prefix"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])

