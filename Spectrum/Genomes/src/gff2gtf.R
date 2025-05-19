library(rtracklayer)
library(GenomicRanges)
#system args to get file from nextflow
args = commandArgs(trailingOnly=T)
gff = args[1]
anno <- rtracklayer::import(gff)
#gff <- renameSeqlevels(gff, gff[1])  
rtracklayer::export(anno,"convert2.gtf","gtf")

