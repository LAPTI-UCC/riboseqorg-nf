// workflow for fetching data

include { FASTQ_DL } from '../../modules/local/fastq_dl.nf'
include { FASTQC } from '../../modules/local/fastqc.nf'

workflow fetch_data {

    take: samples_ch

    main:
        fastq_path_ch   =   FASTQ_DL( samples_ch )
        fastqc_path_ch  =   FASTQC( fastq_path_ch )

    emit:
        fastq_path_ch
}
