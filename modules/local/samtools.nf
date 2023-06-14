
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
	file sorted_bam

	output:
	path "${sorted_bam.baseName}.bam_sorted", emit: sorted_bam
	path "${sorted_bam.baseName}.bam_sorted.bai", emit: sorted_bam_bai
	

	"""
	samtools index ${sorted_bam.baseName}.bam_sorted
	"""

}