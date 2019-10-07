import "split_by_region.wdl" as SBR #workflow for the inside of the file splitting loop
import "create_indexes.wdl" as CI #workflow for the index generation before merging

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
	String python_dir

        command {
                python ${python_dir}get_ch_len.py ${fasta_file} ${Ch}
        }
        output {
                String ch_len = read_string(stdout())
        }
}

task get_regions {
        Int max
        Int overlap
        Int region
	String python_dir

        command {
                python ${python_dir}get_regions.py ${max} ${overlap} ${region}
        }
        output {
                Array[Array[String]] regions = read_tsv(stdout())
        }
}

task select_files_by_region {
        Array[String] region
	Array[File] files
	String python_dir	

        command {
                python ${python_dir}select_files_by_region.py  "${region[0]}-${region[1]}" ${sep=" " files} 
        }

        output {
               Array[Array[File]] files_by_region = read_tsv(stdout())
        }
}

task generate_region_file {
	Array[String] reg
	Array[File] files
	Array[File] indexes
	Array[String] bnames

	command {
		mv ${sep=" " files} .
		bcftools index ${sep="; bcftools index " bnames}
		bcftools merge -o spike.X.${reg[0]}-${reg[1]}.vcf.gz ${sep=" " bnames}
	}

#	output {
#		File abr = "spike.X.${reg[0]}-${reg[1]}.vcf.gz"
#	}
}

workflow OcAsPhaser {
	Array[String] Samples
	String Ch
	File vcf
	String Name
	String bam_dir
        String fasta_file
	String python_dir
        Int region_size
        Int overlap_size

        call get_chromosome_length {input: fasta_file=fasta_file, Ch=Ch, python_dir=python_dir}	
	call get_regions {input: max=get_chromosome_length.ch_len, overlap=overlap_size, region=region_size, python_dir=python_dir}
  
	scatter (s in Samples) {
		call get_recode_vcf {input: Ch=Ch, Ind=s, og_vcf=vcf, Name=Name, working_dir="~/phasing_tests/files"}
		call tweek_sample {input: Sample=s}
		call apply_whatshap {input: Name=Name, Ch=Ch, Sample=tweek_sample.return, bam_dir=bam_dir, recode_vcf=get_recode_vcf.recode_vcf}	
		call SBR.split_by_region {input: regions=get_regions.regions, vcf=apply_whatshap.phased_vcf, ch=Ch, Name=Name, Sample=s}
	}

	scatter (r in get_regions.regions) {		
		call select_files_by_region {input: region=r, files=flatten(split_by_region.regs), python_dir=python_dir}
		call CI.create_indexes {input: vcfs=select_files_by_region.files_by_region[0]}
		call generate_region_file {input: reg=r, files=select_files_by_region.files_by_region[0], indexes = create_indexes.indexes, bnames=create_indexes.bnames}
	
	
#	output {
#		Array[Array[File]] by_sample_by_region = SBR.split_by_region.by_region
#	}	
}
