
process SAMTOOLS_NAME_SORT {

	errorStrategy { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

	input:
	file transcriptome_alignments /// from transcriptome_sams ///

	output:
	file "${transcriptome_alignments.baseName}.bam_sorted" /// into sorted_bams ///

	"""
	samtools sort -m 1G -n -@ 8 ${transcriptome_alignments} > ${transcriptome_alignments.baseName}.bam_sorted
	"""
}


process SAMTOOLS_COORD_SORT {

	errorStrategy { task.attempt <= maxRetries  ? 'retry' :  'ignore' }


	input:
	file transcriptome_alignments

	output:
	file "${transcriptome_alignments.baseName}.bam_sorted" 

	"""
	samtools sort -m 1G -@ 8 ${transcriptome_alignments} > ${transcriptome_alignments.baseName}.bam_sorted
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