

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


process BOWTIE_TRANSCRIPTOME {

	publishDir "$params.study_dir/trips_alignment_stats", mode: 'copy', pattern: '*_trips_alignment_stats.txt' 
	publishDir "$params.study_dir/trips_alignments", mode: 'copy', pattern: '*.bam_sorted' 

	input:    
	file less_rrna_fastq /// from fastq_less_rRNA ///

	output:
	path "${less_rrna_fastq.baseName}_transcriptome.bam_sorted", emit: transcriptome_bam
	path "${less_rrna_fastq.baseName}_trips_alignment_stats.txt", emit: mRNA_alignment_stats

	"""
	bowtie -p 8 --norc -a -m 100 -l 25 -n 2 $params.transcriptome_index -q ${less_rrna_fastq} -S 2>> ${less_rrna_fastq.baseName}_trips_alignment_stats.txt |

	samtools view -@ 8 -b -S  | 

	samtools sort -m 1G -@ 8 -o ${less_rrna_fastq.baseName}_transcriptome.bam_sorted 2> blank 1> blank_2
	"""
} 


process BOWTIE_GENOME {

	publishDir "$params.study_dir/gwips_alignment_stats", mode: 'copy', pattern: '*_gwips_alignment_stats.txt'
	publishDir "$params.study_dir/gwips_alignments", mode: 'copy', pattern: '*.bam_sorted'
	
    input:
   	file less_rrna_fastq /// from fastq_less_rRNA ///

    output:
    path "${less_rrna_fastq.baseName}_genome.bam_sorted", emit: genome_sorted_bam /// into genome_sams ///
    path "${less_rrna_fastq.baseName}_gwips_alignment_stats.txt", emit: gwips_alignment_stats  /// into gwips_alignment_stats ///

    """
	bowtie -p 8 -m 1 -n 2 --seedlen 25 ${params.genome_index} -q ${less_rrna_fastq} -S 2>> ${less_rrna_fastq.baseName}_gwips_alignment_stats.txt | 

	samtools view -@ 8 -b -S  |

	samtools sort -m 1G -@ 8 -o ${less_rrna_fastq.baseName}_genome.bam_sorted 2>&1
	"""
}