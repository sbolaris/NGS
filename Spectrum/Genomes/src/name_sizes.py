#take in two files the first is names and second is sizes
import pandas as pd
import sys

names = pd.read_csv(sys.argv[-2], sep="\t", header=None)
sizes = pd.read_csv(sys.argv[-1],  sep="\t", header=None)

merge = pd.concat([names,sizes], axis=1) 

merge.to_csv("genome_sizes", sep="\t", index=False, header=None)
