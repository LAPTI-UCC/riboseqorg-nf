

// Include the necessary modules for this workflow
include { FASTQC } from '../../modules/local/fastqc.nf'

workflow quality_control {
    take:
        fastq

    main:
        fastqc_ch   =   FASTQC(fastq)
    
    emit:
        fastqc_ch
}