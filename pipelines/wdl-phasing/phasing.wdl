import "OcAsPhaser.wdl" as OAP #Workflow for the inside of the whatshap loop

workflow phasing {
	Array[String] Samples
	String og_vcf	
	Array[String] Chrs
	String Name
	String bam_dir
	String fasta_file
	Int region_size
	Int overlap_size
	String File_dir

	# pre-phase all the files through whatshap
	scatter (ch in Chrs){
		call OAP.OcAsPhaser {input: Ch=ch, Samples=Samples, vcf=og_vcf, Name=Name, bam_dir=bam_dir, fasta_file=fasta_file, region_size=region_size, overlap_size=overlap_size, File_dir=File_dir}		
	}

	# obtain each chromosome max_site from fasta file

	# split them on smaller regions with overlap
	# fuse all the files by region
	
	# phase all the files through shapeit4

	# fuse all the files into one big file again
}
