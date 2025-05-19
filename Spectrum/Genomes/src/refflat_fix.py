import pandas as pd 
import sys

ref = sys.argv[1]

tab = pd.read_csv(ref, sep="\t", header=None)

names = tab[0].copy()
tab.insert(loc=1,column="name", value=names)
tab.to_csv("updated_annotations.refflat",header=None, index=False, sep="\t")


