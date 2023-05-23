// workflow for fetching data

include { FASTQ_DL } from '../../modules/local/fastq_dl.nf'

workflow fetch_data {

    take: samples_ch

    main:
        fastq_path_ch   =   FASTQ_DL( samples_ch )

    emit:
        fastq_path_ch
}
