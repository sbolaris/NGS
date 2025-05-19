import pandas as pd
import sys
import re
import csv
import ray

gtf = sys.argv[1]
#cores = sys.argv[2]

ray.init()

tab = pd.read_csv(gtf, sep="\t", header=None)

@ray.remote
def rename_gene(index):
    try:                                                                     
        gene_name = re.search('gene_name "(.+?)";', tab.iloc[i,8]).group(1)  
        tab.iloc[i,0] = gene_name                                            
    except AttributeError:                                                   
        pass                                                                 

result_ids = []

for i in range(0, len(tab[8])):
    result_ids.append(rename_gene.remote(i))

results = ray.get(result_ids) 

tab.to_csv("updated_anno.gtf", sep="\t", index=False,header=None,quoting = csv.QUOTE_NONE)

