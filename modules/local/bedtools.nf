process BEDTOOLS_GENOMECOV {

	input:
	file genome_sorted_bam /// genome_aligned_and_sorted_bam ///

	output:
	path "${genome_sorted_bam.baseName}.sorted.cov", emit: sorted_beds

	"""
	bedtools genomecov -ibam ${genome_sorted_bam.baseName}.bam_sorted -g $params.chrom_sizes_file -bg > ${genome_sorted_bam.baseName}.cov
	sort -k1,1 -k2,2n ${genome_sorted_bam.baseName}.cov > ${genome_sorted_bam.baseName}.sorted.cov
	"""
}