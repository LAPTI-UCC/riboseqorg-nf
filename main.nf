#!/usr/bin/env nextflow

/// Specify to use Nextflow DSL version 2
nextflow.enable.dsl=2

/// Import modules and subworkflows
include { quality_control } from './subworkflows/local/quality_control.nf'
inlcude { fetch_data } from './subworkflows/local/fetch_data.nf'

// Log the parameters
log.info """\

=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||                 RiboSeqOrg Data Processing Pipeline                            
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||  Parameters                                                             
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||  Sample Sheet    : ${params.sample_sheet}}                                     
||  outDir          : ${params.output_dir}                                        
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

/// Define the main workflow
workflow {
    /// Define the input channels
    fastq_ch = Channel.fromPath("${params.input_dir}/*.fastq.gz")
                        .ifEmpty { exit 1, "No fastq files found in ${params.input_dir}" }

    /// Run the subworkflow
    fetch_data(params.sample_sheet)
}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}
