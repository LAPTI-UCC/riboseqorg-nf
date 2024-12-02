process GENOME_BAM_TO_BED {

    input:
	file genome_index_sorted_bam /// from genome_bams ///
	file genome_index_sorted_bam_bai
	val offset

    output:
	path "${genome_index_sorted_bam.baseName}.bam_sorted.sorted.bed", emit: sorted_beds /// into sorted_beds ///
    	
    """
	python3 $projectDir/scripts/bam_to_bed.py ${genome_index_sorted_bam.baseName}.bam_sorted $offset  $params.genome_fasta
	sort -k1,1 -k2,2n ${genome_index_sorted_bam.baseName}.bam_sorted.bed > ${genome_index_sorted_bam.baseName}.bam_sorted.sorted.bed
	"""
}