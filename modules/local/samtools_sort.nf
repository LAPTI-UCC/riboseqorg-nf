
process SAMTOOLS_NAME_SORT {

	input:
	file transcriptome_sam /// from transcriptome_sams ///

	output:
	file "${transcriptome_sam.baseName}.bam_sorted" /// into sorted_bams ///

	"""
	samtools view -@ 8 -b -S ${transcriptome_sam.baseName}.sam -o ${transcriptome_sam.baseName}.bam
	samtools sort -m 1G -n -@ 8 ${transcriptome_sam.baseName}.bam > ${transcriptome_sam.baseName}.bam_sorted
	"""
}