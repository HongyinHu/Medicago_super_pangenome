# 使用interproscan进行注释

import re
import os 
import sys

script = "path/to/software/interproscan/interproscan-5.30-69.0/interproscan.sh"

def main(pep_fa):
	cmd = f"{script} -i {pep_fa} -cpu 20 -iprlookup -goterms -f TSV -dp -pa 2>>run.log"
	print(cmd)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s pep_fa"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])