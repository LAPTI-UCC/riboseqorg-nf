#!/usr/bin/env nextflow

params.sra_files = "./sra/*.sra"

sra_files_ch = Channel.fromPath(params.sra_files)

process sra_to_fastq {

	input:
	file sra_files from sra_files_ch

	output:
	file '*.fastq' into fastq_files_channel1,fastq_files_channel2

	"""
	fastq-dump $sra_files 
	"""
}


process fastqc_on_raw {
	publishDir 'fastqc_on_raw', mode: 'copy'
	
	input:
	file raw_fastq from fastq_files_channel1

	output:
	file '*_fastqc.{zip,html}' into raw_fastqc_dir

	"""
	fastqc -q $raw_fastq 
	"""
}


process clip_fastq {
        
        input:
        file raw_fastq from fastq_files_channel2

        output:
        file '*_clipped.fastq' into clipped_fastq_channel
	
	script:
	if( params.adapter1 && params.adapter2)
        	"""
        	cutadapt --minimum-length=25 -a $params.adapter1 -a $params.adapter2 -o $raw_fastq"_clipped.fastq" $raw_fastq
        	"""
	else
       		"""
		cutadapt --minimum-length=25 -a $params.adapter1 -o $raw_fastq"_clipped.fastq" $raw_fastq
		"""
}


process rRNA_mapping {
	publishDir 'less_rRNA_fastq_files', mode: 'copy', pattern: '*_less_rRNA.fastq'
	publishDir 'rRNA_alignment_stats', mode: 'copy', pattern: '*_rRNA_stats.txt'

	input: 
	file clipped_fastq from clipped_fastq_channel

	output:
	file "${clipped_fastq.baseName}_rRNA_stats.txt" into rRNA_stats
	file "${clipped_fastq.baseName}_less_rRNA.fastq" into fastq_less_rRNA_ch1,fastq_less_rRNA_ch2

	"""
	bowtie -p 8 -v 3 --norc --phred33-qual $params.rRNA_index -q ${clipped_fastq} --un ${clipped_fastq.baseName}_less_rRNA.fastq > ${clipped_fastq.baseName}_rRNA_stats.txt 2>&1
	"""
}




if (params.skip_trips == false){
	process transcriptome_mapping {
		publishDir 'trips_alignment_stats', mode: 'copy', pattern: '*_trips_alignment_stats.txt' 

		input:    
		file less_rrna_fastq from fastq_less_rRNA_ch1

		output:
		file "${less_rrna_fastq.baseName}_transcriptome.sam" into transcriptome_sams
		file "${less_rrna_fastq.baseName}_trips_alignment_stats.txt" into mRNA_alignment_stats

		"""
		bowtie -p 8 --norc -a -m 100 -l 25 -n 2  -S  -x $params.transcriptome_index -q ${less_rrna_fastq} ${less_rrna_fastq.baseName}_transcriptome.sam  > ${less_rrna_fastq.baseName}_trips_alignment_stats.txt 2>&1
		"""
	} 


	process transcriptome_sam_to_bam {
		input:
		file transcriptome_sam from transcriptome_sams

		output:
		file "${transcriptome_sam.baseName}.bam_sorted" into sorted_bams

		"""
		samtools view -@ 8 -b -S ${transcriptome_sam.baseName}.sam -o ${transcriptome_sam.baseName}.bam
		samtools sort -m 1G -n -@ 8 ${transcriptome_sam.baseName}.bam > ${transcriptome_sam.baseName}.bam_sorted
		"""
	}




	process bam_to_sqlite {
		publishDir 'sqlites', mode: 'copy', pattern: '*.sqlite'
		input:
		file sorted_bam from sorted_bams

		output:
		file "*.sqlite" into sqlite_ch

		"""
		bam_to_sqlite.py ${sorted_bam} $params.annotation_sqlite ${sorted_bam.baseName}
		"""
	}

}



if (params.skip_gwips == false) {
	process genome_mapping {
		publishDir 'gwips_alignment_stats', mode: 'copy', pattern: '*_gwips_alignment_stats.txt'

	        input:
        	file less_rrna_fastq from fastq_less_rRNA_ch2

	        output:
        	file "${less_rrna_fastq.baseName}_genome.sam" into genome_sams
        	file "${less_rrna_fastq.baseName}_gwips_alignment_stats.txt" into gwips_alignment_stats

	        """
        	bowtie -p 8 -m 1 -n 2 --seedlen 25 -S -x ${params.genome_index} -q ${less_rrna_fastq} -S ${less_rrna_fastq.baseName}_genome.sam  >> ${less_rrna_fastq.baseName}_gwips_alignment_stats.txt 2>&1
        	"""
	}




	process genome_sam_to_bed {
        	input:
        	file genome_sam from genome_sams

	        output:
        	file "${genome_sam.baseName}.sorted.cov" into coverage_beds
		file "${genome_sam.baseName}.bam_sorted.sorted.bed" into sorted_beds
       		"""
        	samtools view -@ 8 -b -S ${genome_sam.baseName}.sam -o ${genome_sam.baseName}.bam
        	samtools sort -m 1G -@ 8 ${genome_sam.baseName}.bam > ${genome_sam.baseName}.bam_sorted
		samtools index ${genome_sam.baseName}.bam_sorted
		bam_to_bed.py ${genome_sam.baseName}.bam_sorted 15  $params.genome_fasta
		sort -k1,1 -k2,2n ${genome_sam.baseName}.bam_sorted.bed > ${genome_sam.baseName}.bam_sorted.sorted.bed
		bedtools genomecov -ibam ${genome_sam.baseName}.bam_sorted -g $params.chrom_sizes_file -bg > ${genome_sam.baseName}.cov
		sort -k1,1 -k2,2n ${genome_sam.baseName}.cov > ${genome_sam.baseName}.sorted.cov
	        """
	}



	process bed_to_bigwig {
		publishDir 'bigwigs', mode: 'copy', pattern: '*.bw'

		input:
		file bedfile from sorted_beds

		output:
		file "*.bw" into bigwigs

		"""
		bedGraphToBigWig ${bedfile} $params.chrom_sizes_file ${bedfile.baseName}.bw
		"""
	}


	process coveragebed_to_bigwig {
        	publishDir 'bigwigs', mode: 'copy', pattern: '*.bw'

	        input:
        	file bedfile from coverage_beds

	        output:
        	file "*.bw" into cov_bigwigs

        	"""
        	bedGraphToBigWig ${bedfile} $params.chrom_sizes_file ${bedfile.baseName}.coverage.bw
        	"""
	}
}






