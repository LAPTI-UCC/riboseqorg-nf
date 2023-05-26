// workflow for preprocessing data 

include { FIND_ADAPTERS } from '../../modules/local/find_adapters.nf'
include { CUTADAPT } from '../../modules/local/cutadapt.nf'

workflow preprocessing {

    take: 
        fastq_ch
        samples_ch
    

    main:
        adapter_ch          =   FIND_ADAPTERS   ( fastq_ch )
        trimmed_fastq_ch    =   CUTADAPT        ( fastq_ch, adapter_ch)

    emit:
        trimmed_fastq_ch
}
