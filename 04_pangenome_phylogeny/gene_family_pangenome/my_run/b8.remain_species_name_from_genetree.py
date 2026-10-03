# 删除基因树中的基因，保留物种名称

import re
import os 
import sys

def dels(i):
	res = re.search(r"(|.*?):",i).group(1)
	return res

def main(file):
	with open(file) as inf:
		for lines in inf.readlines():
			line = re.sub(r"\|\w+[.]*[-]*\w+[.]*\w+[.]*\w+","",lines)
			print(line,end="")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("python %s all_tree_file"%sys.argv[0])
        sys.exit(0)
    main(sys.argv[1])

