using GenomicAnnotations
chr = readgbk("/hdd-roo/DHS/gbk/Pseudomonas_syringae_MAFF212063_6757.gbk")[1]

for gene in chr.genes
	gene.locus_tag = "$(chr.name)_$(gene.locus_tag)"
	println(gene.locus_tag)
end
