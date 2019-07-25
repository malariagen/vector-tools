import "split_by_region.wdl" as SBR #workflow for the inside of the file splitting loop


task get_recode_vcf {
        String Ch
        String Ind
        File og_vcf
        String Name
	String working_dir

        command {
                vcftools --gzvcf ${og_vcf} --chr ${Ch} --indv ${Ind} --out ${Name}.${Ch}.${Ind} --recode
        }
        output {   
                File recode_vcf = "${Name}.${Ch}.${Ind}.recode.vcf"
        }
}

task apply_whatshap {
	String Name
	String Ch
	String Sample
	String bam_dir
	File recode_vcf

	command {
		whatshap phase -o ${Name}.${Ch}.${Sample}.phased.vcf ${recode_vcf} ${bam_dir}${Sample}.bam
	}
	
	output {
		File phased_vcf = "${Name}.${Ch}.${Sample}.phased.vcf"
	}
}

task tweek_sample {
        String Sample

        command {
		echo ${Sample} | tr - _
        }
        output {
                String return = read_string(stdout())
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

workflow OcAsPhaser {
	Array[String] Samples
	String Ch
	File vcf
	String Name
	String bam_dir
        String fasta_file
        Int region_size
        Int overlap_size


	scatter (s in Samples) {
		call get_recode_vcf {input: Ch=Ch, Ind=s, og_vcf=vcf, Name=Name, working_dir="~/phasing_tests/files"}
		call tweek_sample {input: Sample=s}
#		call apply_whatshap {input: Name=Name, Ch=Ch, Sample=tweek_sample.return, bam_dir=bam_dir, recode_vcf=get_recode_vcf.recode_vcf} 
                call get_chromosome_length {input: fasta_file=fasta_file, Ch=Ch}
                call get_regions {input: max=get_chromosome_length.ch_len, overlap=overlap_size, region=region_size}		
#		call SBR.split_by_region {input: regions=get_regions.regions, vcf=apply_whatshap.phased_vcf, ch=Ch, Name=Name, Sample=s}
		call SBR.split_by_region {input: regions=get_regions.regions, vcf=get_recode_vcf.recode_vcf, ch=Ch, Name=Name, Sample=s}
	}
	
#	output {
#		Array[File] recode_vcfs = apply_whatshap.phased_vcf
#	}	
}
