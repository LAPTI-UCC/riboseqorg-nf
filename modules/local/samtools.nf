
process SAMTOOLS_NAME_SORT {

	errorStrategy { task.attempt <= maxRetries  ? 'retry' :  'ignore' }


	input:
	file transcriptome_alignments /// from transcriptome_sams ///

	output:
	file "${transcriptome_alignments.baseName}.bam_sorted" /// into sorted_bams ///

	"""
	samtools sort -m 1G -n -@ 8 ${transcriptome_alignments} > ${transcriptome_alignments.baseName}.bam_sorted
	"""
}

process SAMTOOLS_INDEX {

	errorStrategy { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

	input:
	file genome_sorted_bam

	output:
	path "${genome_sorted_bam.baseName}.bam_sorted", emit: genome_index_sorted_bam ///not outputting the index///
	path "${genome_sorted_bam.baseName}.bam_sorted.bai", emit: genome_index_sorted_bam_bai
	

	"""
	samtools index ${genome_sorted_bam.baseName}.bam_sorted
	"""

}