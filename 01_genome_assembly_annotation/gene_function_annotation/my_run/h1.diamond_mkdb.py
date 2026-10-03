# 使用diamond进行比对，数据库建库
# 数据库可以是压缩文件

import re
import os 
import sys

script = "path/to/home/anaconda3/envs/biosofeware/bin/diamond"

def main(database):
	name = os.path.basename(database)
	cmd = f"{script} makedb --in {database} -d {name}"
	print(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s database_pep"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])