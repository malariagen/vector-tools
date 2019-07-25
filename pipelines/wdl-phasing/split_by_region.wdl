task get_region {
	File vcf
	String ch
	Int init
	Int stop
	String Name
	String Sample

	command {
		bgzip -c ${vcf} > ${Name}.${ch}.${Sample}.gz
		tabix -fp vcf ${Name}.${ch}.${Sample}.gz
		tabix -h ${Name}.${ch}.${Sample}.gz ${ch}:${init}-${stop} > ${Name}.${ch}.${Sample}.${init}-${stop}.vcf 
	}
	output {
		File vcf_r = "${Name}.${ch}.${Sample}.${init}-${stop}.vcf"
	}
}

workflow split_by_region {
	Array[Array[String]] regions
	File vcf
	String ch
        String Name
        String Sample

	scatter (r in regions) {
		call get_region {input: vcf=vcf, ch=ch, init=r[0], stop=r[1], Name=Name, Sample=Sample}
	}
}
