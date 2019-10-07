import "OcAsPhaser.wdl" as OAP #Workflow for the inside of the whatshap loop

task get_samples_bam { #copies the bam files
	String Sample
	
	command{
		cp /kwiat/vector/mirror/bam/bwa_gatk/${Sample}.bam .
	}
	output {
		File out = "${Sample}.bam"
	}
}

task get_samples_bai { #copies the bai files
	String Sample

	command {
		cp /kwiat/vector/mirror/bam/bwa_gatk/${Sample}.bam.bai .
	}
	output {
		File out = "${Sample}.bam.bai"
	}
}

task get_og_vcf { #copies the original vcf file
	File og_vcf

	command {
		cp ${og_vcf} .
	}
	output {
		File vcf = "${og_vcf}"
		File tbi = "${og_vcf}.tbi"
	}
}

task gen_list_samples_indv { #generate the list of indvs for vcftools/whatshap --unused
	Array[String] Samples

	command {
		Lsi=""
		for s in ${sep= ' ' Samples} ; do
			Lsi=$Lsi" --indv $s"
		done
		echo $Lsi >> /home/jon/phasing_tests/echo.txt
	}
	output {
		String out = "$Lsi" 
	}
}

task get_genome_region { #uses tabix to extract only part of the vcf files
	File og_vcf
	File og_vcf_tbi
	String Region
	String Name

	command {
		tabix -h ${og_vcf} ${Region} > ${Name}.vcf
	}
	output {
		File genome_region = "${Name}.vcf"
	}
}

task get_recode_vcf { #uses vcftools to recode the vcf file before sending it to whatshap 
        String Ch
        String Ind
        File og_vcf
        String Name

        command {
                vcftools --gzvcf ${og_vcf} --chr ${Ch} --indv ${Ind} --out ${Name}.${Ch}.${Ind} --recode
        }
        output {   
                File recode_vcf = "${Name}.${Ch}.${Ind}.recode.vcf"
        }
}

task tweek_sample { #replace the '-' in sample names with a '_'
	String Sample

	command {
		echo ${Sample} | tr - _
	}
	output {
		String return = read_string(stdout())
	}
}

task get_fasta_file { #copies the fasta file that will be used to import chromosome lengths
	String fasta_file

	command {
		cp ${fasta_file} fasta.fa
	}
	output {
		File fa = "fasta.fa"
	}
}

task get_chromosome_length { #from the fasta file, extract the length of each chromosome
	String Ch
	String fasta_file

	command {
		python ~/phasing_tests/get_ch_len.py ${fasta_file} ${Ch}
	} 
	output {
		String ch_len = read_string(stdout())
	}
}

task get_regions {
	Int max
	Int overlap
	Int region

	command {
		python ~/phasing_tests/get_regions.py ${max} ${overlap} ${region}
	}
	output {
		Array[Array[String]] regions = read_tsv(stdout())
	}
}



workflow phasing {
	Array[String] Samples
	String og_vcf	
	Array[String] Chrs
	String Name
	String bam_dir
	String fasta_file
	Int region_size
	Int overlap_size
	String python_dir

	#get original vcf
	#call tweek_sample {input: Sample="AC0091-C"}

	#get the fasta file
	#call get_fasta_file {input: fasta_file=fasta_file} 

	# split it by chromosome and sample
	# pre-phase all the files through whatshap
	scatter (ch in Chrs){
		call OAP.OcAsPhaser {input: Ch=ch, Samples=Samples, vcf=og_vcf, Name=Name, bam_dir=bam_dir, fasta_file=fasta_file, python_dir=python_dir, region_size=region_size, overlap_size=overlap_size}		
	}

	# obtain each chromosome max_site from fasta file

	# split them on smaller regions with overlap
	# fuse all the files by region
	
	# phase all the files through shapeit4

	# fuse all the files into one big file again
}
