# 获取两两物种间共线性
# 参考物种在前，自己物种放后
# 突出显示某一条，simple 文件中前面加上g*， g为green

import re
import os
import sys

def main(species1, species2):
	print("python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names {a} {b}".format(a=species1,b=species2))
	print("python -m jcvi.compara.synteny screen --minspan=30 --simple {a}.{b}.anchors {a}.{b}.anchors.new".format(a=species1,b=species2))
	# print("python -m jcvi.graphics.karyotype {seqids} {layout}".format(seqids=seqids, layout=layout))

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s species1 species2"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1], sys.argv[2])
