task get_region {
	File vcf
	String ch
	Int init
	Int stop
	String Name
	String Sample
	String file_dir


	command {
		tabix -h ${file_dir}${Name}.${ch}.${Sample}.gz ${ch}:${init}-${stop} > ${file_dir}${Name}.${ch}.${Sample}.${init}-${stop}.vcf
		bgzip -c ${file_dir}${Name}.${ch}.${Sample}.${init}-${stop}.vcf > ${file_dir}${Name}.${ch}.${Sample}.${init}-${stop}.vcf.gz
		rm ${file_dir}${Name}.${ch}.${Sample}.${init}-${stop}.vcf
	}
	output {
		File vcf_r = "${file_dir}${Name}.${ch}.${Sample}.${init}-${stop}.vcf.gz"
	}
}

workflow split_by_region {
	Array[Array[String]] regions
	File vcf
	String ch
        String Name
        String Sample
	String file_dir

	scatter (r in regions) {
		call get_region {input: vcf=vcf, ch=ch, init=r[0], stop=r[1], Name=Name, Sample=Sample, file_dir=file_dir}
	}
}
