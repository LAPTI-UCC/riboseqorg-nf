process BOWTIE_RRNA {

	publishDir "$params.study_dir/less_rRNA_fastq_files", mode: 'copy', pattern: '*_less_rRNA.fastq'
	publishDir "$params.study_dir/rRNA_alignment_stats", mode: 'copy', pattern: '*_rRNA_stats.txt'

	input: 
	file clipped_fastq /// from clipped_fastq_channel ///

	output:
	path "${clipped_fastq.baseName}_rRNA_stats.txt" , emit: rRNA_stats
	path "${clipped_fastq.baseName}_less_rRNA.fastq", emit: fastq_less_rRNA

	"""
	bowtie -p 8 -v 3 --norc --phred33-qual $params.rRNA_index -q ${clipped_fastq} --un ${clipped_fastq.baseName}_less_rRNA.fastq 2> ${clipped_fastq.baseName}_rRNA_stats.txt 
	"""
}