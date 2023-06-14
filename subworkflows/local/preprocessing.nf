// workflow for preprocessing data 

include { FIND_ADAPTERS } from '../../modules/local/find_adapters.nf'
include { CUTADAPT } from '../../modules/local/cutadapt.nf'
include { COLLAPSE_FASTQ } from '../../modules/local/collapse_fastq.nf'
include { BOWTIE_RRNA } from '../../modules/local/bowtie.nf'

workflow preprocessing {

    take: 
        fastq_ch
        samples_ch
    

    main:
        adapter_ch          =   FIND_ADAPTERS   ( fastq_ch )
        trimmed_fastq_ch    =   CUTADAPT        ( fastq_ch, adapter_ch)
        collapsed_fastq_ch  =   COLLAPSE_FASTQ  ( trimmed_fastq_ch )

    emit:
        collapsed_fastq_ch
}
