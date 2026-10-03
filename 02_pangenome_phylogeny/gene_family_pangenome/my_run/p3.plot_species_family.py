# 绘制每个物种中的分类分布

import re
import os 
import sys
import pandas as pd
from plotnine import *

def main(file):
	data = pd.read_csv(file,sep="\t",header=0)
	df = pd.melt(data,id_vars="Sample")
	category = ["Core","SoftCore","Shell","Specific"]
	df["variable"] = df["variable"].astype("category").cat.set_categories(category)
	Sort_df = df.sort_values(by="value",ascending=False)
	print(Sort_df)

	p1 = ggplot(Sort_df,aes(x="Sample",y="value",fill='variable')) + \
		geom_bar(stat="identity",color="black",position="fill",width=0.6,size=0.25) + \
		scale_fill_manual(values={"Core":"#E36A6E","SoftCore":"#E58DA2","Shell":"#5F9DD4","Specific":"#E6DA90"}) + \
		coord_flip() + \
		theme_classic()

	p1.save("p3.barplot_specific_family.pdf",dpi=600,width=6,height=6)

main(sys.argv[1])