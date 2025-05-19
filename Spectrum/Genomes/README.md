# How to make a genome for Sequoia Projects

## Required files 
1. Genome - NCBI genomes page right click and copy link for genome in fasta format  wget on the link that says 
2. GFF / GTF file - NCBI genomes page right click copy link address for the gff link in the upperleft of the page followed by a wget to save the file to local machine 
3. Annotation Report - NCBI genome page click on link for for assembly under related information, wget the link "Full sequence report"

## Create genome\_path
1. save above files to a folder of what ever you would like to call it
1. unzip files with gunzip (typically gff/gtf and the genome file (.fna)
- if the genome ends with .fas you will need to change the extension to .fna otherwise your pipeline will think its done with out processing the genome  
1. point the nextflow option of --genome\_path  to that folder

## Full Run Command:
```
nextflow run ~/repos/genome-annotations/nextflow/make\_genome\_index.nf --genome\_name [your name for the genome] --genome\_path [dir with genome and annoation files]
```

# Easy Way

## Ensemble 

If your genome is on ensemble http://uswest.ensembl.org/info/data/ftp/index.html then congrats you can use the genome and GTF here and it will be high quality and you wont have to fix bad gff files from NCBI or other databases. 

### format changes 
The genome you find on ensemble will match the GTF file all you need to do is wget on the link for the genome and after uncompressing change file extenstion of .fna 
For the GTF it is in the file tree under genes directory 
  
# programatic fixes (TODO)

## Last min fixes
There were some last min fixes to this to make the genomes added and preprocessing to update 

### ref files
Long gtf file - needs a key word gene "gene name here"; this is done with the long\_gtf\_fix.py 
refflat file - the gene pred program that produces the reffflat file is missing a column which needs to have the name of the gene this is fixed with the refflat\_fix.py
biotypes - This is fixed now but it needed the gene name that matches the key word gene in the long gtf file - however this produces no biotype for ercc (see below) 

### ERCC addition 
There is an issue with nextflow that if you concat to a variable that symlink is updated through out the pipeline to avoid this making a file once and adding the ercc after wards is preferable.
spike gtf: this can be added to the contig gtf
spike long gtf: this can be added to a long gtf (has the key word gene already)
biotypes: already set up to have the biotype of ercc

