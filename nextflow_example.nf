#!/usr/bin/env nextflow


println "Script dir: $baseDir"
println "Working dir: $workDir"

params.raw_reads = "<absolute_path_to_ribo_dir>/*_rpf_1.fastq.gz"
raw_fastq_datasets = Channel.fromPath(params.raw_reads)
raw_fastq_datasets.into {datasets_fastqc_on_raw; datasets_cutadapt}

process fastqc_on_raw {
	publishDir '<absolute_path_to_ribo_dir>/fastqc_on_raw', mode: 'copy'
	
	input:
	file raw_reads from datasets_fastqc_on_raw

	output:
	file 'fastqc_on_raw/*_fastqc.{zip,html}' into raw_fastqc_dir

	"""
	mkdir -p <absolute_path_to_ribo_dir>/fastqc_on_raw
	mkdir -p fastqc_on_raw
	fastqc -q $raw_reads -o fastqc_on_raw
	"""
}

/*
 * Trim adapter sequence with cutadapt
 */
 

adapter_1 = 'AAAAAAAA'

process cutadapt_remove_adapter {
	publishDir '<absolute_path_to_ribo_dir>/clip_adapter', mode: 'copy'
	
	input:
	file raw_reads from datasets_cutadapt

	output:
	file "${raw_reads.baseName}" into clip_adapter_fastq_dir

	"""		
	mkdir -p <absolute_path_to_ribo_dir>/clip_adapter
	cutadapt -a $adapter_1 --minimum-length=23 ${raw_reads} -o ${raw_reads.baseName} --untrimmed-output <absolute_path_to_ribo_dir>/clip_adapter/${raw_reads.baseName}_untrimmed --too-short-output <absolute_path_to_ribo_dir>/clip_adapter/${raw_reads.baseName}_tooshort > <absolute_path_to_ribo_dir>/clip_adapter/${raw_reads.baseName}_cutadapt_log.txt
	"""
}

clip_adapter_fastq_dir.into {datasets_fastx_trimmer}

process fastx_trim {
	publishDir '<absolute_path_to_ribo_dir>/trimmed_reads', mode: 'copy'
	
	input:
	file clipped_reads from datasets_fastx_trimmer

	output:
	file "${clipped_reads.baseName}" into trimmed_fastq_dir

	"""
	mkdir -p <absolute_path_to_ribo_dir>/trimmed_reads
	mkdir -p trimmed_reads
	fastx_trimmer -Q33 -f 4 -i ${clipped_reads} -o ${clipped_reads.baseName}
	"""
}


trimmed_fastq_dir.into {datasets_bowtie_rRNA}

/*
 * Map rRNA reads using bowtie
 */
process rRNA_mapping {
	publishDir '<absolute_path_to_ribo_dir>/mapped_rRNA', mode: 'copy'
	
	input: 
	file clipped_fastq from datasets_bowtie_rRNA

	output:
	file "${clipped_fastq.baseName}_rRNA.sam" into mapped_rRNA_dir
	file "${clipped_fastq.baseName}_rRNA_stats.txt" into rRNA_stats
	file "${clipped_fastq.baseName}_less_rRNA" into fastq_less_rRNA_dir

	"""
	mkdir -p <absolute_path_to_ribo_dir>/mapped_rRNA
	mkdir -p <absolute_path_to_ribo_dir>/less_rRNA
	bowtie -v 3 --norc -S --phred33-qual -q <absolute_path_to_rRNA_index> ${clipped_fastq} ${clipped_fastq.baseName}_rRNA.sam --un ${clipped_fastq.baseName}_less_rRNA > ${clipped_fastq.baseName}_rRNA_stats.txt 2>&1
	cp ${clipped_fastq.baseName}_less_rRNA <absolute_path_to_ribo_dir>/less_rRNA
	"""
}

fastq_less_rRNA_dir.into {datasets_bowtie_tRNA}

/*
 * Map tRNA reads using bowtie
 */
process tRNA_mapping {
	publishDir '<absolute_path_to_ribo_dir>/mapped_tRNA', mode: 'copy'
	
	input: 
	file clipped_fastq from datasets_bowtie_tRNA

	output:
	file "${clipped_fastq.baseName}_tRNA.sam" into mapped_tRNA_dir
	file "${clipped_fastq.baseName}_tRNA_stats.txt" into tRNA_stats
	file "${clipped_fastq.baseName}_less_tRNA" into fastq_less_tRNA_dir

	"""
	mkdir -p <absolute_path_to_ribo_dir>/mapped_tRNA
	mkdir -p <absolute_path_to_ribo_dir>/less_tRNA
	bowtie -v 3 --norc -S --phred33-qual -q <absolute_path_to_tRNA_index> ${clipped_fastq} ${clipped_fastq.baseName}_tRNA.sam --un ${clipped_fastq.baseName}_less_tRNA > ${clipped_fastq.baseName}_tRNA_stats.txt 2>&1
	cp ${clipped_fastq.baseName}_less_tRNA <absolute_path_to_ribo_dir>/less_tRNA
	"""
}



fastq_less_tRNA_dir.into {datasets_fastqc_less_tRNA; datasets_bowtie_mRNA}

process fastqc_less_tRNA {
	publishDir '<absolute_path_to_ribo_dir>/fastqc_less_tRNA', mode: 'copy'
	
	input:
	file raw_reads from datasets_fastqc_less_tRNA

	output:
	file 'fastqc_less_tRNA/*_fastqc.{zip,html}' into fastqc_less_tRNA_dir

	"""
	mkdir -p <absolute_path_to_ribo_dir>/fastqc_less_tRNA
	mkdir -p fastqc_less_tRNA
	fastqc -q $raw_reads -o fastqc_less_tRNA
	"""
}


/*
 * Map reads using bowtie - to the transcriptome in this case for 
 */
process mRNA_mapping { 
	publishDir '<absolute_path_to_ribo_dir>/Transcriptome_alignments/alignments', mode: 'copy'
	   
	input:    
	file fastq_less_tRNA from datasets_bowtie_mRNA

	output:
	file "${fastq_less_tRNA.baseName}.sam" into mapped_mRNA_dir
	file "${fastq_less_tRNA.baseName}_mRNA_alignment_stats.txt" into mRNA_alignment_stats
	//file "${fastq_less_tRNA.baseName}_unmapped" into unmapped_dir

	"""
	mkdir -p <absolute_path_to_ribo_dir>/Transcriptome_alignments/alignments
	mkdir -p <absolute_path_to_ribo_dir>/Transcriptome_alignments/Unmapped
	bowtie -a -m 100 -l 25 -n 2 --norc -S --phred33-qual -q <absolute_path_to_bowtie_index> ${fastq_less_tRNA} ${fastq_less_tRNA.baseName}.sam --un <absolute_path_to_ribo_dir>/Transcriptome_alignments/Unmapped/${fastq_less_tRNA.baseName}_unmapped > ${fastq_less_tRNA.baseName}_mRNA_alignment_stats.txt 2>&1
	"""
} 

mapped_mRNA_dir.into {datasets_for_sam_accepted23; datasets_for_bam}

params.gene_details = "<absolute_path_to_gene_details>"
Human_gene_details = Channel.fromPath(params.gene_details)

/*
 * separate reads that map once, twice, trice 
 */
process sam_accepted23 {
	input:
	file sam_file from datasets_for_sam_accepted23 
	file gene_details from Human_gene_details.first()
	
	output:
	file "${sam_file.baseName}_1.sam" into accepted_alignments1
	file "${sam_file.baseName}_2.sam" into accepted_alignments2 
	file "${sam_file.baseName}_3.sam" into accepted_alignments3 

	script:
	"""
	mkdir -p <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments1
	mkdir -p <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments2
	mkdir -p <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments3
	Patrick_sam_accepted23 ${gene_details} ${sam_file}
	cp ${sam_file.baseName}_1.sam <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments1
	cp ${sam_file.baseName}_2.sam <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments2
	cp ${sam_file.baseName}_3.sam <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments3
	"""
}


/*
 * Convert bowtie output in sam format to bam format
 */
process Convert_sam_to_bam {
	input:
	file original_sam_file from datasets_for_bam
	file accepted_sam_file1 from accepted_alignments1
	file accepted_sam_file2 from accepted_alignments2
	file accepted_sam_file3 from accepted_alignments3

	output:
	file "${original_sam_file.baseName}.bam" into bam_dir
	file "${accepted_sam_file1.baseName}.bam" into bam_dir1
	file "${accepted_sam_file2.baseName}.bam" into bam_dir2
	file "${accepted_sam_file3.baseName}.bam" into bam_dir3

	script:
	"""
	samtools view -Sb ${original_sam_file} | samtools sort > ${original_sam_file.baseName}.bam 
	samtools view -Sb ${accepted_sam_file1} | samtools sort > ${accepted_sam_file1.baseName}.bam 
	samtools view -Sb ${accepted_sam_file2} | samtools sort > ${accepted_sam_file2.baseName}.bam 
	samtools view -Sb ${accepted_sam_file3} | samtools sort > ${accepted_sam_file3.baseName}.bam
	cp ${original_sam_file.baseName}.bam <absolute_path_to_ribo_dir>/Transcriptome_alignments/alignments
	cp ${accepted_sam_file1.baseName}.bam <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments1
	cp ${accepted_sam_file2.baseName}.bam <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments2
	cp ${accepted_sam_file3.baseName}.bam <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments3
	"""
}

process index_bam {

	input:
	file bam_file from bam_dir
	file bam_file1 from bam_dir1
	file bam_file2 from bam_dir2
	file bam_file3 from bam_dir3

	output:
	file "${bam_file.baseName}.bam.bai" into bam_dir_index
	file "${bam_file1.baseName}.bam.bai" into bam_dir_index1
	file "${bam_file2.baseName}.bam.bai" into bam_dir_index2
	file "${bam_file3.baseName}.bam.bai" into bam_dir_index3

	script:
	"""
	samtools index ${bam_file}
	samtools index ${bam_file1}
	samtools index ${bam_file2}
	samtools index ${bam_file3}
	cp ${bam_file.baseName}.bam.bai <absolute_path_to_ribo_dir>/Transcriptome_alignments/alignments
	cp ${bam_file1.baseName}.bam.bai <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments1
	cp ${bam_file2.baseName}.bam.bai <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments2
	cp ${bam_file3.baseName}.bam.bai <absolute_path_to_ribo_dir>/Transcriptome_alignments/accepted_alignments3
	"""
}


