// workflow for preprocessing data 

include { FIND_ADAPTERS } from '../../modules/local/find_adapters.nf'
include { FASTP } from '../../modules/local/fastp.nf'
include { FASTQC } from '../../modules/local/fastqc.nf'
include { COLLAPSE_FASTQ } from '../../modules/local/collapse_fastq.nf'

workflow preprocessing {

    take: 
        fastq_ch
        samples_ch
    

    main:
        fastqc_ch           =   FASTQC          ( fastq_ch )
        adapter_ch          =   FIND_ADAPTERS   ( fastq_ch, fastqc_ch.fastqc_data )

        trimmed_fastq_ch    =   FASTP           ( fastq_ch, adapter_ch )

        collapsed_fastq_ch  =   COLLAPSE_FASTQ  ( trimmed_fastq_ch.trimmed_fastq)

    emit:
        collapsed_fastq_ch
}
