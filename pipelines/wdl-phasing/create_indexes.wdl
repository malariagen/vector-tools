task create_index {
	File vcf
	String bname

	command {
		bcftools index ${vcf} > ${bname}.csi
	}
	output {
		File ind = "${bname}.csi"
	}
}

workflow create_indexes {
	Array[File] vcfs

	scatter (f in vcfs) {
		String tbname = basename(f)
		call create_index {input: vcf=f, bname=tbname}
	}

	output {
		Array[File] indexes = create_index.ind
		Array[String] bnames = tbname
	}
}
