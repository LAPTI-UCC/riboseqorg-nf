// workflow for fetching data

include { FASTQ_DL } from '../../modules/local/fastq_dl.nf'
include { FASTQC } from '../../modules/local/fastqc.nf'
include { FETCH_RUN } from '../../modules/local/fetch.nf'
include { FASTERQ_DUMP } from '../../modules/local/fasterq_dump.nf'

workflow fetch_data {

    take: samples_ch

    main:
        // fastq_ch   =   FASTQ_DL( samples_ch )
        sra_ch          =   FETCH_RUN       ( samples_ch )
        fastq_ch        =   FASTERQ_DUMP    ( sra_ch )

    emit:
        fastq_ch
        samples_ch
}
