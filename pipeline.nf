/* THE pipeline */

 /* OLD LINE OF CODE params.sra_files = "./sra/*.sra"  */

/* -------------------
PRE-PROCESSING BRANCH
--------------------- */

/*VERY provisional -> I am literally copy and pasting the full path to the file. Ideally, this should be done from terminal when calling nextflow*/

params.fastq_files = ./data/2019_Homo_sapiens_GSE125114_SRP179636/fastq/*.fastq

process clip_fastq {
        
    input:
    file raw_fastq 

    output:
    file '*_clipped.fastq' /* into clipped_fastq_channel  */
	
	script: 
	/* There was an if statement here, referring to two adapter paratmeters. Since we are using a file now, I deleted that line of code.*/
	"""
    cutadapt --minimum-length=25 -a "file:adapters.fa;min_overlap=5;noindels" -o $raw_fastq"_clipped.fastq" $raw_fastq
    """
}

process rRNA_mapping {
	publishDir 'less_rRNA_fastq_files', mode: 'copy', pattern: '*_less_rRNA.fastq'
	publishDir 'rRNA_alignment_stats', mode: 'copy', pattern: '*_rRNA_stats.txt'

	input: 
	file clipped_fastq /* from clipped_fastq_channel */

	output:
	path "${clipped_fastq.baseName}_rRNA_stats.txt" , emit: rRNA_stats
	path "${clipped_fastq.baseName}_less_rRNA.fastq", emit: fastq_less_rRNA

	"""
	bowtie -p 8 -v 3 --norc --phred33-qual $params.rRNA_index -q ${clipped_fastq} --un ${clipped_fastq.baseName}_less_rRNA.fastq > ${clipped_fastq.baseName}_rRNA_stats.txt 2>&1
	"""
}

/* ORIGINALLY THE BELOW PROCESS WAS NAMED "fastqc_on_raw". It has been updated for consistency, considering we are
using fastqc on processed reads in this new version (sequences with no adapters and no rRNAs)-> new name is fastqc_on_processed */

process fastqc_on_processed {
	publishDir 'fastqc_on_processed', mode: 'copy'
	
	input:
	file processed_fastq /*should I call it fastq_less_rRNA instrad? doesn't change a thing technically, but yh,know*/

	output:
	file '*_fastqc.{zip,html}' /* into raw_fastqc_dir */

	"""
	fastqc -q $processed_fastq
	"""
}

process multiqc_on_fastq {

	input:
	/* should specify the directory? above there is a publishDir thing */
	file ('fastqc/*')

	output:
	file "multiqc_report.html"
	file "multiqc_output_data"

	/* script to execute fastqc*/
	"""
	multiqc .
	"""
}


/* -------------------------
TRANSCRIPTOME MAPPING BRANCH
---------------------------- */


process transcriptome_mapping {
	publishDir 'trips_alignment_stats', mode: 'copy', pattern: '*_trips_alignment_stats.txt' 

	input:    
	file less_rrna_fastq /* from fastq_less_rRNA */

	output:
	path "${less_rrna_fastq.baseName}_transcriptome.sam", emit: transcriptome_sams  /* USE AN EMIT COMMAND HERE?*/
	path "${less_rrna_fastq.baseName}_trips_alignment_stats.txt", emit: mRNA_alignment_stats

	"""
	bowtie -p 8 --norc -a -m 100 -l 25 -n 2  -S  -x $params.transcriptome_index -q ${less_rrna_fastq} ${less_rrna_fastq.baseName}_transcriptome.sam  > ${less_rrna_fastq.baseName}_trips_alignment_stats.txt 2>&1
	"""
} 

process transcriptome_sam_to_bam {
	input:
	file transcriptome_sam /* from transcriptome_sams */

	output:
	file "${transcriptome_sam.baseName}.bam_sorted" /* into sorted_bams */

	"""
	samtools view -@ 8 -b -S ${transcriptome_sam.baseName}.sam -o ${transcriptome_sam.baseName}.bam
	samtools sort -m 1G -n -@ 8 ${transcriptome_sam.baseName}.bam > ${transcriptome_sam.baseName}.bam_sorted
	"""
}

process bam_to_sqlite {
	publishDir 'sqlites', mode: 'copy', pattern: '*.sqlite'
	input:
	file sorted_bam /* from sorted_bams */

	output:
	file "*.sqlite" /* into sqlite_ch */

	"""
	bam_to_sqlite.py ${sorted_bam} $params.annotation_sqlite ${sorted_bam.baseName}
	"""
}


/* --------------------
GENOME MAPPING BRANCH 
----------------------*/


	
process genome_mapping {
	publishDir 'gwips_alignment_stats', mode: 'copy', pattern: '*_gwips_alignment_stats.txt'
    input:
   	file less_rrna_fastq /* from fastq_less_rRNA */

    output:
    path "${less_rrna_fastq.baseName}_genome.sam", emit: genome_sams /* into genome_sams */
    path "${less_rrna_fastq.baseName}_gwips_alignment_stats.txt", emit: gwips_alignment_stats/* into gwips_alignment_stats */

    """
	bowtie -p 8 -m 1 -n 2 --seedlen 25 -S -x ${params.genome_index} -q ${less_rrna_fastq} -S ${less_rrna_fastq.baseName}_genome.sam  >> ${less_rrna_fastq.baseName}_gwips_alignment_stats.txt 2>&1
	"""
}


process genome_sam_to_bed {
    input:
	file genome_sam /* from genome_sams */

    output:
    path "${genome_sam.baseName}.sorted.cov", emit: coverage_beds /* into coverage_beds */
	path "${genome_sam.baseName}.bam_sorted.sorted.bed", emit: sorted_beds /* into sorted_beds */
    	
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
	file bedfile /* from sorted_beds */
	
    output:
	file "*.bw"  /* into bigwigs */

	"""
	bedGraphToBigWig ${bedfile} $params.chrom_sizes_file ${bedfile.baseName}.bw
	"""
}


process coveragebed_to_bigwig {
	publishDir 'bigwigs', mode: 'copy', pattern: '*.bw'

	input:
    file bedfile /* from coverage_beds */

	output:
    file "*.bw"  /* into cov_bigwigs */

	"""
	bedGraphToBigWig ${bedfile} $params.chrom_sizes_file ${bedfile.baseName}.coverage.bw
	"""
}







/* THE WORKFLOW BLOCK: It specifies the order of the processes and where outputs are used as inputs*/

workflow {
    fastq_data = Channel.fromPath(params.fastq_files) /* Assign the fastq files in a folder to fastq_data */
    clip_fastq(fastq_data)   /* Uses fastq_data and clips away the adapters */
	rRNA_mapping(clip_fastq.out)
	fastqc_on_processed(rRNA_mapping.out.fastq_less_rRNA)
    multiqc_on_fastq(fastqc_on_processed.out)		

    /* IF STATEMENT #1 */
    if (params.skip_trips == false) {
        transcriptome_mapping(rRNA_mapping.out.fastq_less_rRNA)
        transcriptome_sam_to_bam(transcriptome_mapping.out.transcriptome_sams)
        bam_to_sqlite(transcriptome_sam_to_bam.out)
    }
    /* IF STATEMENT #2 */
    if (params.skip_gwips == false) {
        genome_mapping(rRNA_mapping.out.fastq_less_rRNA)
        genome_sam_to_bed(genome_mapping.out.genome_sams)
        bed_to_bigwig(genome_sam_to_bed.out.sorted_beds)
        coveragebed_to_bigwig(genome_sam_to_bed.out.coverage_beds)
    }
}
