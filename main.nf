#!/usr/bin/env nextflow

/// Specify to use Nextflow DSL version 2
nextflow.enable.dsl=2

/// Import modules and subworkflows
include { fetch_data } from './subworkflows/local/fetch_data.nf'
include { STAR_ALIGN } from './modules/local/STAR/main.nf'
include { RIBOMETRIC } from './modules/local/ribometric/main.nf'
include { SAMTOOLS_COORD_SORT } from './modules/local/samtools/samtools_coord_sort/main.nf'
include { SAMTOOLS_NAME_SORT } from './modules/local/samtools/samtools_name_sort/main.nf'
include { BAM_TO_SQLITE } from './modules/local/bam_to_sqlite/main.nf'

// include { BOWTIE_RRNA } from './modules/local/bowtie.nf'

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
||  outdir       : ${params.outdir}                                     
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
    samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            def meta = [
                id: row.Run,
                study_accession: row.study_accession,
            ]
            [ meta, row.Run ]
        }
    fetch_data_ch           =   fetch_data(samples_ch)

    STAR_ALIGN(
        fetch_data_ch,
        params.star_index,
        params.gtf
    )
    SAMTOOLS_COORD_SORT(
        STAR_ALIGN.out.transcriptome_bam,
    )
    RIBOMETRIC(
        SAMTOOLS_COORD_SORT.out.bam,
        params.ribometric_annotation
    )

    SAMTOOLS_NAME_SORT(STAR_ALIGN.out.transcriptome_bam)
    BAM_TO_SQLITE(
        SAMTOOLS_NAME_SORT.out.bam,
        params.annotation_sqlite
    )

}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}
