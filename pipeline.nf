/* THE pipeline */

/* OLD LINE OF CODE params.sra_files = "./sra/*.sra" */

/* -------------------
PRE-PROCESSING BRANCH
--------------------- */

/*VERY provisional -> I am literally copy and pasting the full path to the file. Ideally, this should be done from terminal when calling nextflow*/

project_dir = projectDir  /*specify a new variable, the project directory */

process CLIP_FASTQ {
        
    input:
    file raw_fastq 

    output:
    file '*_clipped.fastq' /* into clipped_fastq_channel  */
	
	script: 
	"""
    cutadapt --minimum-length=25 -a "file:$params.adapter_fasta" -o $raw_fastq"_clipped.fastq" $raw_fastq
    """
}

process rRNA_MAPPING {

	publishDir "$params.study_dir/less_rRNA_fastq_files", mode: 'copy', pattern: '*_less_rRNA.fastq'
	publishDir "$params.study_dir/rRNA_alignment_stats", mode: 'copy', pattern: '*_rRNA_stats.txt'

	input: 
	file clipped_fastq /* from clipped_fastq_channel */

	output:
	path "${clipped_fastq.baseName}_rRNA_stats.txt" , emit: rRNA_stats
	path "${clipped_fastq.baseName}_less_rRNA.fastq", emit: fastq_less_rRNA

	"""
	bowtie -p 8 -v 3 --norc --phred33-qual $params.rRNA_index -q ${clipped_fastq} --un ${clipped_fastq.baseName}_less_rRNA.fastq 2> ${clipped_fastq.baseName}_rRNA_stats.txt 
	"""
}

/* ORIGINALLY THE BELOW PROCESS WAS NAMED "fastqc_on_raw". It has been updated for consistency, considering we are
using fastqc on processed reads in this new version (sequences with no adapters and no rRNAs)-> new name is fastqc_on_processed */

process FASTQC_ON_PROCESSED {

	publishDir "$params.study_dir/fastqc", mode: 'copy'
	
	input:
	file processed_fastq

	output:
	file '*_fastqc.{zip,html}' /* into raw_fastqc_dir */

	"""
	fastqc -q $processed_fastq
	"""
}

process MULTIQC_ON_FASTQ {

	publishDir "$params.study_dir/multiqc", mode: 'copy'

	input:
	file ('fastqc/*')

	output:
	file "multiqc_report.html"
	
	"""
	multiqc $params.study_dir/fastqc
	"""
}


/* -------------------------
TRANSCRIPTOME MAPPING BRANCH
---------------------------- */


process TRANSCRIPTOME_MAPPING {

	publishDir "$params.study_dir/trips_alignment_stats", mode: 'copy', pattern: '*_trips_alignment_stats.txt' 

	input:    
	file less_rrna_fastq /* from fastq_less_rRNA */

	output:
	path "${less_rrna_fastq.baseName}_transcriptome.sam", emit: transcriptome_sam
	path "${less_rrna_fastq.baseName}_trips_alignment_stats.txt", emit: mRNA_alignment_stats

	"""
	bowtie -p 8 --norc -a -m 100 -l 25 -n 2 $params.transcriptome_index -q ${less_rrna_fastq} -S ${less_rrna_fastq.baseName}_transcriptome.sam  > ${less_rrna_fastq.baseName}_trips_alignment_stats.txt 2>&1
	"""
} 

process TRANSCRIPTOME_SAM_TO_BAM {

	input:
	file transcriptome_sam /* from transcriptome_sams */

	output:
	file "${transcriptome_sam.baseName}.bam_sorted" /* into sorted_bams */

	"""
	samtools view -@ 8 -b -S ${transcriptome_sam.baseName}.sam -o ${transcriptome_sam.baseName}.bam
	samtools sort -m 1G -n -@ 8 ${transcriptome_sam.baseName}.bam > ${transcriptome_sam.baseName}.bam_sorted
	"""
}

process BAM_TO_SQLITE {

	publishDir "$params.study_dir/sqlites", mode: 'copy', pattern: '*.sqlite'

	input:
	file sorted_bam /* from sorted_bams */

	output:
	file "*.sqlite" /* into sqlite_ch */

	"""
	python3 $project_dir/scripts/bam_to_sqlite.py ${sorted_bam} $params.annotation_sqlite ${sorted_bam.baseName}
	"""
}


/* --------------------
GENOME MAPPING BRANCH 
----------------------*/

process GENOME_MAPPING {

	publishDir "$params.study_dir/gwips_alignment_stats", mode: 'copy', pattern: '*_gwips_alignment_stats.txt'
	
    input:
   	file less_rrna_fastq /* from fastq_less_rRNA */

    output:
    path "${less_rrna_fastq.baseName}_genome.bam_sorted", emit: genome_sorted_bam /* into genome_sams */
    path "${less_rrna_fastq.baseName}_gwips_alignment_stats.txt", emit: gwips_alignment_stats/* into gwips_alignment_stats */

    """
	bowtie -p 8 -m 1 -n 2 --seedlen 25 ${params.genome_index} -q ${less_rrna_fastq} -S 2>> ${less_rrna_fastq.baseName}_gwips_alignment_stats.txt | 

	samtools view -@ 8 -b -S |

	samtools sort -m 1G -@ 8 -o ${less_rrna_fastq.baseName}_genome.bam_sorted
	"""
}



process INDEX_SORT_BAM {

	input:
	file genome_sorted_bam

	output:
	path "${genome_sorted_bam.baseName}.bam_sorted", emit: genome_index_sorted_bam


	"""
	samtools index ${genome_sorted_bam.baseName}.bam_sorted
	"""

}


process BAM_TO_COVBED {

	input:
	file genome_sorted_bam /// this should be changed to make reading easier ///

	output:
	path "${genome_sorted_bam.baseName}.sorted.cov", emit: coverage_beds

	"""
	bedtools genomecov -ibam ${genome_sorted_bam.baseName}.bam_sorted -g $params.chrom_sizes_file -bg > ${genome_sorted_bam.baseName}.cov
	sort -k1,1 -k2,2n ${genome_sorted_bam.baseName}.cov > ${genome_sorted_bam.baseName}.sorted.cov
	"""
}


process GENOME_BAM_TO_BED {

    input:
	file genome_index_sorted_bam /// from genome_bams ///

    output:
	path "${genome_index_sorted_bam.baseName}.bam_sorted.sorted.bed", emit: sorted_beds /// into sorted_beds ///
    	
    """
	python3 $project_dir/scripts/bam_to_bed.py ${genome_index_sorted_bam.baseName}.bam_sorted 15  $params.genome_fasta
	sort -k1,1 -k2,2n ${genome_index_sorted_bam.baseName}.bam_sorted.bed > ${genome_index_sorted_bam.baseName}.bam_sorted.sorted.bed
	"""
}


process BED_TO_BIGWIG {

	publishDir "$params.study_dir/bigwigs", mode: 'copy', pattern: '*.bw'

	input:
	file bedfile /* from sorted_beds */
	
    output:
	file "*.bw"  /* into bigwigs */

	"""
	$project_dir/scripts/bedGraphToBigWig ${bedfile} $params.chrom_sizes_file ${bedfile.baseName}.bw
	"""
}

/* NEED TO:
1) perform a sanity check over the code (expecially name of I/O)
2) modify the workflow object accordingly... since the condition is still not implemented, maybe just update it and write a draft
3) check how the file name is changed across the various processes. we want the final beds from each branch to have different names 


/*
process COVERAGEBED_TO_BIGWIG {

	publishDir "$params.study_dir/bigwigs", mode: 'copy', pattern: '*.bw'

	input:
    file bedfile /* from coverage_beds */

///	output:
///   file "*.bw"  /* into cov_bigwigs */

///	"""
///	$project_dir/scripts/bedGraphToBigWig ${bedfile} $params.chrom_sizes_file ${bedfile.baseName}.coverage.bw
///	"""

///} 


workflow {

	fastq_data = Channel.fromPath ( params.fastq_files )

	CLIP_FASTQ          ( fastq_data )
	rRNA_MAPPING        ( CLIP_FASTQ.out )
	FASTQC_ON_PROCESSED ( rRNA_MAPPING.out.fastq_less_rRNA )
	MULTIQC_ON_FASTQ    ( FASTQC_ON_PROCESSED.out )		

    /// TRANSCRIPTOME MAPPING ///
	if ( params.skip_trips == false ) {
		
		TRANSCRIPTOME_MAPPING    ( rRNA_MAPPING.out.fastq_less_rRNA )
		TRANSCRIPTOME_SAM_TO_BAM ( TRANSCRIPTOME_MAPPING.out.transcriptome_sams )
		BAM_TO_SQLITE            ( TRANSCRIPTOME_SAM_TO_BAM.out )

	}

    /// GENOME MAPPING ///
	if ( params.skip_gwips == false ) {

		GENOME_MAPPING        ( rRNA_MAPPING.out.fastq_less_rRNA )
		/// This block is for RNA-Seq studies only. It's executed depending on a parameter, which defines the type of study we are working with.
		params.x = "something_temporary"
		if (params.x == "Is a RNA-Seq study") {
			BAM_TO_COVBED     ( GENOME_MAPPING.out.genome_sorted_bam )
			BED_TO_BIGWIG	  ( BAM_TO_COVBED.out.coverage_beds)
		}
		/// The following block is executed if the study is not an RNA-Seq. ///
		else {
			INDEX_SORT_BAM    ( GENOME_MAPPING.out.genome_sorted_bam )
			GENOME_BAM_TO_BED ( INDEX_SORT_BAM.out.genome_index_sorted_bam )
			BED_TO_BIGWIG     ( GENOME_BAM_TO_BED.out.sorted_beds )
			BAM_TO_COVBED     ( INDEX_SORT_BAM.out.genome_index_sorted_bam )
			BED_TO_BIGWIG     ( GENOME_BAM_TO_BED.out.coverage_beds )
		}

		/// OLD CODE ///
		/*
		GENOME_BAM_TO_BED     ( GENOME_MAPPING.out.genome_sorted_bams )
		BED_TO_BIGWIG         ( GENOME_BAM_TO_BED.out.sorted_beds )
		COVERAGEBED_TO_BIGWIG ( GENOME_BAM_TO_BED.out.coverage_beds )
		*/
    }

}

