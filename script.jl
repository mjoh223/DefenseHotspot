using GenomicAnnotations

genbank_files = readdir("/hdd-roo/DHS/gbk/", join=true)
for file in genbank_files
	records = readgbk(file)
	for record in records
		for gene in record.genes
			gene.locus_tag = "$(record.name)_$(gene.locus_tag)"
			println(gene.locus_tag)
		end
	end
end
