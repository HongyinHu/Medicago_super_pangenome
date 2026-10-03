# 重新命名gff文件中的geneid名称

import re
import os 
import sys

pre_spec = "Vvi"
pub_geneid = ""
num = 1

def main(gff):
	with open(gff) as inf:
		for lines in inf.readlines():
			if not lines.startswith("#"):
				line = lines.strip().split()

				if line[2] == "gene":
					pub_geneid = "%s%06d"%(pre_spec,num)
					num +=1
					line[-1] = "ID=%s"%pub_geneid

				if line[2] == ""



