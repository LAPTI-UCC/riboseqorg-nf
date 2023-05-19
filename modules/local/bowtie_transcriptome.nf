// This should pass straight to samtools sort 

process BOWTIE_TRANSCRIPTOME {

	publishDir "$params.study_dir/trips_alignment_stats", mode: 'copy', pattern: '*_trips_alignment_stats.txt' 

	input:    
	file less_rrna_fastq /// from fastq_less_rRNA ///

	output:
	path "${less_rrna_fastq.baseName}_transcriptome.sam", emit: transcriptome_sam
	path "${less_rrna_fastq.baseName}_trips_alignment_stats.txt", emit: mRNA_alignment_stats

	"""
	bowtie -p 8 --norc -a -m 100 -l 25 -n 2 $params.transcriptome_index -q ${less_rrna_fastq} -S ${less_rrna_fastq.baseName}_transcriptome.sam  > ${less_rrna_fastq.baseName}_trips_alignment_stats.txt 2>&1
	"""
} 