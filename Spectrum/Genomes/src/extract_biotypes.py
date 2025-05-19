import pandas as pd
import sys 
import re
import ray

ray.init()

gtf = sys.argv[1]

tab = pd.read_csv(gtf, sep="\t", header=None) 
#need the model_name, Type,
#tab = tab[[0,8]]
#change headers to gene_id and gene_biotype
biotype = pd.DataFrame()
types = []
gene_names = []


@ray.remote
def get_biotype(i,row):
    #types = []
    #gene_names =[]
    try:
        name = re.search('gene_name "(.+?)"',row[8]).group(1)
        t = re.search('gene_biotype "(.+?)"',row[8]).group(1)
        #types = types+ [t]
        #gene_names = gene_names + [name]
    except AttributeError:
        try:
            name = re.search('gene_name "(.+?)"',row[8]).group(1)
            t = re.search('Note "(.+?)"',row[8]).group(1)   
            #types = types +[t] 
            #gene_names = gene_names +[name]
        except AttributeError: 
            try:
                name = re.search('gene_name "(.+?)"',row[8]).group(1)
                t = re.search('transcript_biotype "(.+?)"',row[8]).group(1)
            except AttributeError:
                try:
                    name = re.search('gene_id "(.+?)"',row[8]).group(1)
                    t = "No Biotype"
                except AttributeError:
                    name = "NO GENE"
                    t = "NO TYPE"
            #types = types +["No Biotype"]
            #gene_names = gene_names +[name]

    return(name,t)

result_ids = []
for i in range(0,len(tab[8])):    
    result_ids.append(get_biotype.remote(i,tab.iloc[i]))

results = ray.get(result_ids)
biotype = pd.DataFrame(results)
biotype.columns = ["gene_id", "gene_biotype"]
#biotype["gene_id"] = gene_names
#biotype["gene_biotype"] = types

biotype = biotype.drop_duplicates()

biotype.to_csv("gene_biotypes.tsv", sep="\t", index=False)
