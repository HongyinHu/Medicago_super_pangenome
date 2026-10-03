# 绘制基因家族中物种数量分布条形图

import re
import os 
import sys
import pandas as pd
from plotnine import *
import matplotlib.pyplot as plt

def main(frequency_family,frequency_bing):
	data1 = pd.read_csv(frequency_family,sep=r"\s+",names=["freq_sample","Fam_num","class"])
	data2 = pd.read_csv(frequency_bing,sep=r"\s+",names=["class","Fam_num"])

	p1 = ggplot(data1,aes(x="freq_sample",y="Fam_num",fill="class"))+ geom_bar(stat="identity",width=0.6) + \
	scale_fill_manual(values={"Core":"#E36A6E","Shell":"#5F9DD4","SoftCore":"#E58DA2","Specific":"#E6DA90"}) + theme_classic()
	p1.save("p2.barplot_frequency_family.pdf",dpi=300,width=6,height=6)
	
	color_dict = {"Core":"#E36A6E","Shell":"#5F9DD4","softCore":"#E58DA2","Specific":"#E6DA90"}
	color_list = [color_dict[label] for label in data2["class"]]
	
	plt.figure(figsize=(6,6))
	plt.pie(data2["Fam_num"],labels=data2["class"],autopct='%1.1f%%',colors=color_list)
	plt.savefig("p2.bingplot_frequency_family.pdf")

	
	# p2.save("p2.bingplot_frequency_family.pdf",dpi=300,width=6,height=6)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print("python %s frequency_family frequency_bing_file"%sys.argv[0])
		sys.exit(0)
	main(sys.argv[1],sys.argv[2])