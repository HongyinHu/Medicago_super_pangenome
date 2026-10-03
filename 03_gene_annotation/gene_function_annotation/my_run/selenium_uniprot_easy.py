#根据创建网址，爬取uniprot基因描述信息
#Edit by hhy, 2023.12.6

import re
import os
import sys
import requests
from bs4 import BeautifulSoup

def make_url(geneid="gene1",accession="O76080"):
	#创建网址并抓取html内容
	url = f"https://rest.uniprot.org/uniprotkb/{accession}.xml"
	content=requests.get(url,timeout=10)
	soup = BeautifulSoup(content.text,'lxml-xml')

	#选取部分信息
	lineage = soup.find("lineage").strings
	if "Viridiplantae" in list(lineage):
		try:
			protein_name = soup.find("recommendedName").fullName.text
		except:
			protein_name = "None"
		
		try:
			gene_name = soup.find("name",type="primary").text
		except:
			gene_name = "None"

		try:
			organism_name = soup.find("name",type="scientific").text
		except:
			organism_name = "None"
		
		try:
			disruption_pheno = soup.find("comment",type="disruption phenotype").text
		except:
			disruption_pheno = "None"
		
		try:
			family_name = soup.find("comment",type="similarity").text
		except:
			family_name = "None"
		
		try:
			finction_des = soup.find("comment",type="function").text
		except:
			finction_des = "None"

		try:
			GO_terms = soup.find_all("dbReference",type="GO")
			GO_terms_list = list()
			for i in GO_terms:
				go_id = re.search(r'id="(.*?)"\s+type="GO"',str(i)).group(1)
				go_term = re.search(r'type="term" value="(.*?)"',str(i)).group(1)
				GO_terms_list.append((go_id+"---"+go_term))
		except:
			GO_terms_list = "None"

		#格式化输出
		geneid_dict = dict()
		geneid_dict.setdefault(geneid,{})["gene_name"] = gene_name.strip()
		geneid_dict.setdefault(geneid,{})["protein_name"] = protein_name.strip()
		geneid_dict.setdefault(geneid,{})["organism_name"] = organism_name
		geneid_dict.setdefault(geneid,{})["family_name"] = "family_name:"+family_name.strip()
		geneid_dict.setdefault(geneid,{})["disruption_pheno"] = "disruption_pheno:"+disruption_pheno.strip()
		geneid_dict.setdefault(geneid,{})["finction_des"] = "function_descript:"+finction_des.strip()
		geneid_dict.setdefault(geneid,{})["GO_term"] = "\t".join(GO_terms_list)

		for k,v in geneid_dict.items():
			return (k+"###"+"###".join(v.values()))
	else:
		print(f"{accession} not in plants")
	
def main(file):
	total,number = sum(1 for _ in open(file)),1

	with open(file) as inf, open(f"{file}.anno.txt","w") as ouf:
		for lines in inf.readlines():
			line = lines.strip().split()
			geneid = line[0]
			accession = line[1].split("|")[1]
			print(f"{number}/{total}/{os.path.basename(file)}--{geneid}/{accession} 正在抓取信息...")
			
			try:
				res = make_url("###".join([line[0],line[1],line[2]]),accession)
				if res:
					print(res,file=ouf)
			except:
				print("###".join([line[0],line[1],line[2],"not get information"]),file=ouf)
			
			number += 1
			

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s uniport_format6_res"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1])