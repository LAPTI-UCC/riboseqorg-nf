#!/usr/bin/env nextflow

/// Specify to use Nextflow DSL version 2
nextflow.enable.dsl=2

include { DATA_ACQUISITION } from './subworkflows/local/data_acquisition'
include { QUALITY_CONTROL } from './subworkflows/local/quality_control'
include { ALIGNMENT } from './subworkflows/local/alignment'
include { POST_PROCESSING } from './subworkflows/local/post_processing'
include { ANALYSIS } from './subworkflows/local/analysis'

// Log the parameters
log.info """\

=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||                 RiboSeqOrg Data Processing Pipeline                            
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||  Parameters                                                             
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||  Sample Sheet    : ${params.sample_sheet}                                     
||  outDir          : ${params.outdir}                                        
||  workDir         : ${workflow.workDir}                             
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

"""
// Help Message to prompt users to specify required parameters
def help() {
    log.info"""
  Usage:  nextflow run main.nf --input <path_to_fastq_dir> 

  Required Arguments:

  --input    Path to directory containing fastq files.

  Optional Arguments:

  --outDir	Path to output directory. 
	
""".stripIndent()
}


workflow {
    // Data Acquisition
    DATA_ACQUISITION(params.sample_sheet)

    // Quality Control
    QUALITY_CONTROL(DATA_ACQUISITION.out.samples)

    // Alignment (only for samples that passed QC)
    ALIGNMENT(
        QUALITY_CONTROL.out.clean_samples,
        params.star_index,
        params.bowtie_index,
        params.rRNA_index,
        params.gtf
        )

    // Post-processing
    POST_PROCESSING(
        ALIGNMENT.out.transcriptome_bam,
        ALIGNMENT.out.genome_bam,
        params.annotation_sqlite,
        params.chrom_sizes_file
    )

    // Analysis
    ANALYSIS(ALIGNMENT.out.transcriptome_bam, params.ribometric_annotation)
}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}
