from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd 
import sys

#get command line arg for file name 
spike_file = sys.argv[-1]
#read in as table 
tab = pd.read_csv(spike_file, sep="\t")
#Loop through spikes to create fasta file
record_list = []
for index, row in tab.iterrows():
    #generate records - id is first tried to keep generic to avoid just ERCC
    record = SeqRecord(Seq(row.Sequence), id=row[0], description="") 
    record_list.append(record)

SeqIO.write(record_list, "spike.fasta", "fasta")
