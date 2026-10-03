# 绘制不同泛基因组家族中kaks的比较箱线图

import re
import os 
import sys
import pandas as pd
from plotnine import *

def main(all_kaks):
	data = pd.read_csv(all_kaks,header=None,sep="\t",names=["a","b"])
	# print(data["a"])
	category = ["Core","SoftCore","Shell","Specific"]
	data.a = data.a.astype("category").cat.set_categories(category)
	# data[[0]] = data[[0]].apply(lambda x: x.astype("category").cat.set_categories(category))
	# data = data.sort_values(by=data[[1]],ascending=False)
	
	p = ggplot(data,aes(x=data.a,y=data.b,fill=data.a)) + \
	geom_boxplot() + \
	scale_fill_manual(values={"Core":"#E36A6E","SoftCore":"#E58DA2","Shell":"#5F9DD4","Specific":"#E6DA90"}) + \
	theme_classic() + ylim(0,2.5)

	p.save("p5.kaks.boxplot.pdf",dpi=600,width=6,height=6)

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print("python %s all_kaks(two columns)")
		sys.exit(0)
	main(sys.argv[1])