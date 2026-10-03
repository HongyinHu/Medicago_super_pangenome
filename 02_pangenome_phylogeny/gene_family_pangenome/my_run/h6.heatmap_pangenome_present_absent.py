# 转换格式绘制家族存在缺失热图

import re
import os
import sys

def main(ortho_count):
	with open(ortho_count) as inf:
		for index,lines in enumerate(inf.readlines()):
			if index > 0:
				line = lines.strip().split()
				for index_l,i in enumerate(line[1:]):
					print(index_l)

main(sys.argv[1])

