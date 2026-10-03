#根据创建网址，爬取NCBI NR基因描述信息
#Edit by hhy, 2023.12.20

import re
import os
import sys
import requests
from bs4 import BeautifulSoup

def make_url(geneid,accession,identity):
	#创建网址并抓取html内容
	url = f"https://www.ncbi.nlm.nih.gov/protein/{accession}"
	content=requests.get(url,timeout=10)
	soup = BeautifulSoup(content.text,'lxml')

	#选取部分信息
	descrip_info = soup.find("div",class_="rprtheader").h1.text.strip()
	return "###".join([geneid,accession,identity,descrip_info])
	
def main(file):
	total,number = sum(1 for _ in open(file)),1

	with open(file) as inf, open(f"{file}.anno.txt","w") as ouf:
		for lines in inf.readlines():
			line = lines.strip().split()
			geneid = line[0]
			accession = line[1]
			print(f"{number}/{total}/{os.path.basename(file)}--{geneid}/{accession} 正在抓取信息...")
			try:
				res = make_url(geneid,accession,line[2])
				print(res,file=ouf)
			except:
				print(line[0],line[1],file=ouf)
			finally:
				number += 1
				

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s uniport_format6_res"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])