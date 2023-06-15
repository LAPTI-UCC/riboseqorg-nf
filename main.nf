#!/usr/bin/env nextflow

/// Specify to use Nextflow DSL version 2
nextflow.enable.dsl=2

/// Import modules and subworkflows
include { quality_control } from './subworkflows/local/quality_control.nf'
include { fetch_data } from './subworkflows/local/fetch_data.nf'
include { preprocessing } from './subworkflows/local/preprocessing.nf'
include { trips_RiboSeq } from './subworkflows/local/trips.nf'
include { gwips_RiboSeq } from './subworkflows/local/gwips.nf'

include { BOWTIE_RRNA } from './modules/local/bowtie.nf'

// Log the parameters
log.info """\

=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||                 RiboSeqOrg Data Processing Pipeline                            
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||  Parameters                                                             
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||  Sample Sheet    : ${params.sample_sheet}                                     
||  outDir          : ${params.output_dir}                                        
||  workDir         : ${workflow.workDir}   
||  study_dir       : ${params.study_dir}                                     
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
    samples_ch  =   Channel
                        .fromPath(params.sample_sheet)
                        .splitCsv(header: true, sep: '\t')
                        .map { row -> tuple("${row.study_accession}", "${row.Run}", "${row.ScientificName}", "${row.LIBRARYTYPE}")}

    fetch_data_ch           =   fetch_data(samples_ch)
    fastq_ch                =   fetch_data_ch.fastq_ch
    samples_ch              =   fetch_data_ch.samples_ch
    collapsed_fastq_ch      =   preprocessing(fastq_ch, samples_ch)
    less_rRNA_ch          =   BOWTIE_RRNA     ( collapsed_fastq_ch )
    if ( params.skip_gwips == false ) {
        gwips_RiboSeq(less_rRNA_ch.fastq_less_rRNA)
    }
    if ( params.skip_trips == false ) {
        trips_RiboSeq(less_rRNA_ch.fastq_less_rRNA)
    }
}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}
