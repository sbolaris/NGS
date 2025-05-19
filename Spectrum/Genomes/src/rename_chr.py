import os
import sys

#args = trailing args
genome = sys.argv[1]
gtf = sys.argv[2]

#grep uniqe chr id from genome 
chr_in_gen = os.popen("grep '>' "+genome+"").read().split('\n')[:-1]
chr_in_gen.sort()
#grep uniqe chr id from gtf 
chr_in_gtf = os.popen("awk -F '\t' '{print $1}' "+gtf+"|uniq").read().split('\n')[:-1] 
#map them 1:1 some how 
for i in range(0,len(chr_in_gen)):
    #use sed to replace the names in genome 
    tmp = chr_in_gen[i]
    key = tmp.split(' ')[0]
    cmd = "sed -i 's/"+key+".*/>"+chr_in_gtf[i]+"/g' "+genome
    os.popen(cmd).read()

