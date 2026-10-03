# 比对结果中筛选每个query的最佳subject
# 只要是format 6格式都可以

import re
import os
import sys

script = "python -m jcvi.formats.blast"

def main(diamond_res):
	cmd = f"conda run -n biosofeware {script} best -n 1 {diamond_res}"
	print(cmd)

if __name__ == '__main__':
	if len(sys.argv) !=2:
		print("python %s diamond_res"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])