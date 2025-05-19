import os
import sys
import pandas as pd
import ray
import re
import csv

ray.init()

long_gtf = sys.argv[1]

#os.popen("grep -v '^#' "+long_gtf+" >tmp.gtf")
tab = pd.read_csv(long_gtf, sep="\t", header=None)

@ray.remote
def add_gene_meta(index, row):
    try:
        name = re.search('gene "(.+?)"',row[8]).group(1)
    except AttributeError:
        try:
            name = re.search('gene_name "(.+?)"',row[8]).group(1)
            gene = ' gene "'+name+'";'
            row[8] = row[8]+gene
        except AttributeError:
            try:
                name = re.search('gene_id "(.+?)"',row[8]).group(1)
                gene = ' gene "'+name+'";'              
                row[8] = row[8]+gene
            except AttributeError: 
                try:
                    name = re.search('transcipt_id "(.+?)"',row[8]).group(1)
                    gene = ' gene "'+name+'";'            
                    row[8] = row[8]+gene                               
                except AttributeError: 
                   print("Failure condition not found")

    return(row)

result_ids = [] 
for i in range(0,len(tab[8])):                          
        result_ids.append(add_gene_meta.remote(i,tab.iloc[i]))
                                                                
results = ray.get(result_ids)                           
long_rna = pd.DataFrame(results)                         
long_rna.to_csv("updated_long.gtf", sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)     
