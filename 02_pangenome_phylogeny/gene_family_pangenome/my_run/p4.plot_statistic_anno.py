# 绘制interpro注释结果

import re
import os 
import sys
from plotnine import *
import pandas as pd

def main(file):
	data = pd.read_excel(file,header=None,engine="openpyxl")
	category = ["Core","SoftCore","Shell","Specific"]
	data.iloc[:,0] = data.iloc[:,0].astype("category").cat.set_categories(category)

	category_fill = ["unanno","anno"]
	data.iloc[:,1] = data.iloc[:,1].astype("category").cat.set_categories(category_fill)
	
	# data_df = data.sort_values(by=data.columns[2],ascending=True)

	p1 = ggplot(data,aes(x=data.iloc[:,0],y=data.iloc[:,2],fill=data.iloc[:,1])) + geom_bar(stat="identity",position="fill") + \
		scale_fill_manual(values={"anno":"#EF7571","unanno":"#22B9B9"}) + theme_classic()
	p1.save("p4.plot_statistic_anno.pdf",dpi=600,width=6,height=6)

main(sys.argv[1])