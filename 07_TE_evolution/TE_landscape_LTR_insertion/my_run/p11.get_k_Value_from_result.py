#从结果文件提取K值

import re
import os
import sys

def main(K_res_indir):
	for file in os.listdir(K_res_indir):
		if file.endswith("meg"):
			fullfile = os.path.join(K_res_indir,file)
			with open(fullfile) as inf:
				content = inf.read()
				# print(content)
				try:
					K_value = re.search(r"\[2\]\s+(\d+.\d+)",content).group(1)
				except:
					continue
				r2 = 5.1*10**(-9)
				# if float(K_value) < 1:
				ltr_time = (float(K_value)/(2*r2))/1000000
				print("\t".join([fullfile,K_value,str(ltr_time)]))

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s K_res_indir"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])


