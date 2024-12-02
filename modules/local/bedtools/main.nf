process BEDTOOLS_GENOMECOV {
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

	input:
	file genome_sorted_bam /// genome_aligned_and_sorted_bam ///

	output:
	path "${genome_sorted_bam.baseName}.sorted.cov", emit: sorted_beds

	"""
	bedtools genomecov -ibam ${genome_sorted_bam.baseName}.bam_sorted -g $params.chrom_sizes_file -bg > ${genome_sorted_bam.baseName}.cov
	sort -k1,1 -k2,2n ${genome_sorted_bam.baseName}.cov > ${genome_sorted_bam.baseName}.sorted.cov
	"""
}