import os
import sys 
import pandas as pd
#import ray

#ray.init()

gtf = sys.argv[-2]
sizes = pd.read_csv(sys.argv[-1],  sep="\t", header=None) 

#@remote
#def get_ribo(index,line):
    #check for Note or gene_biotype as 'rRNA' or pass
    

rint = open('intervals.txt', 'w')
for index, row in sizes.iterrows():
    rint.write("@SQ\tSN:%s\tLN:%s\n" % (row[0], row[1]))

rint.close()

if os.popen("grep 'gene_biotype .rRNA.' "+gtf).read() != '':
    os.popen("grep 'gene_biotype .rRNA.' "+gtf+" | cut -f1,4,5,7,9 >>intervals.txt")

elif os.popen("grep 'transcript_biotype .rRNA.' "+gtf).read() != '':
    os.popen("grep 'transcript_biotype .rRNA.' "+gtf+" | cut -f1,4,5,7,9 >>intervals.txt")

else:
    os.popen("grep 'Note .rRNA.' "+gtf+" | cut -f1,4,5,7,9 >>intervals.txt") 
