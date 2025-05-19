// vim set syntax=groovy

def summary = [:]
summary['Genome Name'] = params.genome_name
summary['Genome Path'] = params.genome_path
summary['GTF file?']   = params.GTFfile
summary['Container']   = params.container
summary['Tar ball']    = params.targz
summary['Cpus']        = params.cpus
summary['RAM']         = params.mem


log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "----------------------------------------------------"


if(params.genome_path == ""){
	println("No Path to genome present")
	println("Please set with --genome_path /path/to/genome/files/")
	help()
	exit 1
}

chr_file = Channel.fromPath("$params.genome_path*.f*a") 
if(params.GTFfile ==false){
	input_anno_file = Channel.fromPath("$params.genome_path*.gff")
}
else{
	input_anno_file = Channel.fromPath("$params.genome_path*.gtf")
}
//spike_file = Channel.fromPath("/opt/biorad/ercc_files/ercc_spike_seq.txt")
//biotype_file = Channel.fromPath("$params.genome_path/*_biotypes.tsv")
anno_report = Channel.fromPath("$params.genome_path/*_assembly_report.txt")
//ercc_gtf = Channel.fromPath("/opt/biorad/ercc_files/ercc_spike.gtf") 
// TODO add preprocessing steps to look for weird stuff such as too much meta data and only take the main info needed for the gtf file from the GFF file 
if(params.GTFfile ==false){
	process covert2gtf{
		publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy'
		
		input:
		file gff_file from input_anno_file
		output:
		file "${params.genome_name}_contigs_in_ref.gtf" into anno_gtf, ribo_gtf, gtf_spike, flat_gtf, longRNA_count_gtf, name_gtf, biotype_gtf
		script:
		//println(gff_file)
		"""
		#need to remove ? from strand and put .
		#check strands of the genes to make sure its all on the same strand
		#check strands shold be + or - if not put .
		sed -i 's/\t?\t/\t.\t/' $gff_file	
		/opt/biorad/src/gffread $gff_file -T -F -o test.gtf
		#python3 /opt/biorad/src/gene_name_update.py test.gtf
		mv test.gtf updated_anno.gtf
		cp updated_anno.gtf ${params.genome_name}_contigs_in_ref.gtf
		"""
	}
}
else{
	process gtf_fn_update{
		publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy'  

		input:	
		file gtf from input_anno_file
		output:
		file "${params.genome_name}_contigs_in_ref.gtf" into anno_gtf,ribo_gtf,gtf_spike,flat_gtf,longRNA_count_gtf,name_gtf,biotype_gtf
		script:
		"""
		cp $gtf ${params.genome_name}_contigs_in_ref.gtf
		""" 
	}
}
/*
process make_biotypes{
	input:
	file gene_anno from biotype_gtf 
	output:
	file ("gene_biotypes.tsv") into biotype_file
	script:
	"""
	grep -v ERCC $gene_anno >real_anno.gtf
	#cp $gene_anno real_anno.gtf
	python3 /opt/biorad/src/extract_biotypes.py real_anno.gtf
	"""

}
*/

//chromosome names of the genome must match gff file 
if(params.GTFfile ==false){  
	process rename_chr{
		publishDir "./$params.genome_name/fasta_ncbi/", mode: 'copy'
		input:
		file chrs from chr_file
		file chr_name from name_gtf
		output:
		file "*.fna" into genome,gen_ercc
		script:
		"""
		python3 /opt/biorad/src/rename_chr.py $chrs $chr_name
		mv $chrs ${params.genome_name}_no_alt_analysis_set.fna 
		"""

	}
}
else{
	process rename_genome{
		publishDir "./$params.genome_name/fasta_ncbi/", mode: 'copy'  
		input:                  
		file chrs from chr_file 
		output:
		file "*.fna" into genome,gen_ercc 
		script:
		"""
		mv $chrs ${params.genome_name}_no_alt_analysis_set.fna
		"""
	}

}
process spike_gtf{
        publishDir "./$params.genome_name/anno_ncbi_ercc/", mode: 'copy'                           
        input:                                                                                     
        file gtf from gtf_spike                                      
        output:                                                                                    
        file("${params.genome_name}_contigs_in_ref.gtf") into anno_gtf_spike, flat_gtf_spike, ribo_gtf_spike, long_gtf_spike       
                                                                                                   
        script:                                                                                    
        """                                                                                        
        cat $gtf >tmp_spike.gtf
	cat tmp_spike.gtf /opt/biorad/ercc_files/ercc_spike.gtf >>${params.genome_name}_contigs_in_ref.gtf                                           
        """                                                                                        	

}
//basic index 
process index{
	publishDir "./$params.genome_name/star_ncbi/", mode: 'copy'
	input:
	file fasta from genome
	file gtf from anno_gtf
	output:
	file "*"
	val "go" into finish_star
	file "chrNameLength.txt" into ribo_sizes
	
	script:
	
	"""
	STAR \
		--runThreadN ${params.cpus} \
		--runMode genomeGenerate \
		--genomeDir ./ \
		--genomeFastaFiles $fasta \
		--sjdbGTFfile $gtf \
		--sjdbOverhang 45
	
	"""
}
process add_spike{          
	publishDir "./$params.genome_name/fasta_ncbi_ercc/", mode: 'copy' 
	input:
	file chr_only from gen_ercc
	//file spike_seq from spike_file
	output:
	file "*_genome_ercc.fna" into genome_ercc, gen_spike_size
	file "spike.fasta" into spike_anno_ercc  
	script:
	"""
	#python script to take txt -> fasta
	python3 /opt/biorad/src/txt2fasta.py /opt/biorad/ercc_files/ercc_spike_seq.txt
	#cat the fastafiles together
	cat $chr_only  >genome_ercc.fna
	cat spike.fasta >>genome_ercc.fna
	mv genome_ercc.fna ${params.genome_name}_genome_ercc.fna
	"""
}
//ercc index
process spike_genome{
	publishDir "./$params.genome_name/star_ncbi_ercc/", mode: 'copy' 
	input:
	file fasta from genome_ercc 
	file gtf from anno_gtf_spike
	val go from finish_star
	output:
	file "*"
	file "chrNameLength.txt" into sizes

	script:
	"""
	STAR \
        --runThreadN ${params.cpus} \
        --runMode genomeGenerate \
        --genomeDir ./ \
        --genomeFastaFiles $fasta \
        --sjdbGTFfile $gtf \
        --sjdbOverhang 45
	"""
}


//sizes 
process get_sizes {         
	publishDir "./$params.genome_name/sizes/", mode: 'copy' 
	input:
	file NameLength from sizes
	
	output:
	file("${params.genome_name}.chrom.sizes")
	file("*_ercc.chrom.sizes") into ribo_sizes_spike

	script:
	"""
	grep -v ^ERCC $NameLength >${params.genome_name}.chrom.sizes
	mv $NameLength ${params.genome_name}_ercc.chrom.sizes 
	"""
}

process gen_ribosomal_int {
	publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy' 
	input:
	file gtf from ribo_gtf
	file sizes from ribo_sizes
	output:
	file("${params.genome_name}_ribosomal_intervals.txt")
	script:
	"""
	python3 /opt/biorad/src/ribosomal_int.py $gtf $sizes
	# might be gene_biotype might be Note
	#grep 'gene_biotype "rRNA"' $gtf | cut -f1,4,5,7,9 >>intervals.txt   
	mv intervals.txt ${params.genome_name}_ribosomal_intervals.txt     
	"""
}

process gen_ribosomal_spike_int {                                                
        publishDir "./$params.genome_name/anno_ncbi_ercc/", mode: 'copy'       
        input:                                                            
        file gtf from ribo_gtf_spike                                            
        file sizes from ribo_sizes_spike                                        
        output:                                                           
        file("${params.genome_name}_ribosomal_intervals.txt")
        script:                                                           
        """                                                               
        python3 /opt/biorad/src/ribosomal_int.py $gtf $sizes
	#grep 'Note "rRNA"' $gtf | cut -f1,4,5,7,9 >>intervals.txt
	mv intervals.txt ${params.genome_name}_ribosomal_intervals.txt 
	"""                                                               
}                                                                         
process gen_flat_file{
	publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy'  
	input:
	file gtf_flat_ref from flat_gtf
	output:
	file("${params.genome_name}_annotations.refflat")
	script:
	"""
	grep -v ERCC $gtf_flat_ref >tmp.gtf
	/opt/biorad/src/gtfToGenePred -allErrors tmp.gtf ${params.genome_name}_annotations.refflat
	python3 /opt/biorad/src/refflat_fix.py ${params.genome_name}_annotations.refflat
	mv updated_annotations.refflat ${params.genome_name}_annotations.refflat 
	"""
}

process gene_flat_spike_file{
	publishDir "./$params.genome_name/anno_ncbi_ercc/", mode: 'copy'
	input:
	file flat_gtf_ref_spike from flat_gtf_spike
	output:
	file("${params.genome_name}_annotations.refflat")
	script: 
	"""                                                                            
	/opt/biorad/src/gtfToGenePred -allErrors $flat_gtf_ref_spike ${params.genome_name}_annotations.refflat   
	python3 /opt/biorad/src/refflat_fix.py ${params.genome_name}_annotations.refflat
	mv updated_annotations.refflat ${params.genome_name}_annotations.refflat
	"""                                                                            
}
process anno_version{
	publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy'
	input:
	file anno_report  
	output:
	file "annotation_version.txt" into anno_v
	script:
	
	"""
	echo "#gtf-version 2.2" >annotation_version.txt
	grep "Assembly name:" $anno_report | awk -F':' '{print \$2}' | xargs echo "#!genome-build " >>annotation_version.txt
	grep "RefSeq assembly accession" $anno_report |awk -F':' '{print \$2}' | xargs echo "#!genome-build-accession " >>annotation_version.txt 
	grep "Date:" $anno_report |awk -F':' '{print \$2}' | xargs echo "#!annotation-date " >>annotation_version.txt  
	grep "Submitter:" $anno_report|awk -F':' '{print \$2}' | xargs echo "#!annotation-source " >>annotation_version.txt 
	"""
}
process anno_version_spike{
	publishDir "./$params.genome_name/anno_ncbi_ercc/", mode: 'copy'
	input:
	file anno_v 
	output:
	file "annotation_version.txt"
	script:
	"""
	echo $anno_v
	"""
}


//need anno_ncbi
process generate_biotypes{
// rename biotype files, generate ref_flat file, and Ribosomal intervals  
	//publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy'
	input:
	file biotype from biotype_gtf
	output:
	file "biotypes.tsv" into spike_biotypes, anno_biotypes
	script:
	"""
	grep -v ^# $biotype >tmp.gtf
	grep -v ^ERCC tmp.gtf >clean.gtf
	python3 /opt/biorad/src/extract_biotypes.py clean.gtf
	mv gene_biotypes.tsv biotypes.tsv
	"""
}
process rename_biotype{
	publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy'
	input:
	file btype from anno_biotypes
	output:
	file "gene_biotypes.tsv"  
	script:
	"""
	cat $btype >gene_biotypes.tsv
	"""
}
//anno_ncbi_ercc
process generate_biotype_spike{
	publishDir "./$params.genome_name/anno_ncbi_ercc/", mode: 'copy'
	input:
	file types from spike_biotypes
	//file ercc_names from spike_anno_ercc //fasta version of the file
	output:
	file "gene_biotypes.tsv"
	script:
	//add ERCC names + ERCC as type 
	"""
	cat $types /opt/biorad/ercc_files/ercc_biotypes.tsv >>gene_biotypes.tsv 
	"""

}

process gen_longRNA_gtf{
	//publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy'
	input:
	file anno_long from longRNA_count_gtf
	output:
	file "updated_long.gtf" into long_gtf_anno, long_gtf_for_spike 
	script:
	"""
	grep -v ^# $anno_long >tmp.gtf
	grep -v ^ERCC tmp.gtf >no_ercc.gtf
	python3 /opt/biorad/src/longRNA_counts_fix.py no_ercc.gtf
	"""
}
process longRNA_gtf{
	publishDir "./$params.genome_name/anno_ncbi/", mode: 'copy' 
	input:
	file long_rna from long_gtf_anno
	output:
	file("${params.genome_name}_longRNA_annotation.gtf")
	script:
	"""
	cat $long_rna >${params.genome_name}_longRNA_annotation.gtf 
	"""
}

process longRNA_spike_gtf{
	publishDir "./$params.genome_name/anno_ncbi_ercc/", mode: 'copy'
	input:
	file long_anno from long_gtf_for_spike 
	//file gtf_anno from ercc_gtf
	output:
	file "${params.genome_name}_longRNA_annotation.gtf"
	script:
	"""
	cat $long_anno >${params.genome_name}_longRNA_annotation.gtf
	cat /opt/biorad/ercc_files/ercc_spike_long.gtf >>${params.genome_name}_longRNA_annotation.gtf  
	"""
}

//create manifests?
/*
process manifest{
	publishDir "./$params.genome_name/", mode: 'copy' 
	input:
	
	ouput:

	script:
	"""
	man md5sum

	"""
}
*/
// TODO : do file collection on all the output files and collect them in a folder to do the tar gz (can take some time)
//tar files for dist if needed
if(params.targz){
	process tarball{
		container 
		publishDir "./$params.genome_name/", mode: 'copy'
		input:
		
		output:
		file "*.tar.gz"
		script:
		"""
		tar -czvf archive_name folder_to_tarball
		"""

	}
} //end if targz

def help(){
	println("nextflow run make_genome_index.nf --genome_name [name] --cpus [cores] --mem [ram] --genome_path [path to files] --GTFfile [true/false]")
	println("The required files from the genome path are the genome with .fna extenstions, gff or gtf file, ans the annotation report of NCBI")
	println("Files can be found by doing NCBI genome search - wget the gnome via link same for gff")
	println("For annoation report, the file is hosted on the assmbly page")

}
