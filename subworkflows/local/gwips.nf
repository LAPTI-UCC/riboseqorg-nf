// workflow for processing files for gwips
include { BOWTIE_GENOME } from '../../modules/local/bowtie.nf'

workflow gwips_RiboSeq {

    take:
        lessRNA_ch    

    main:
        genome_sorted_bam_ch    =   BOWTIE_GENOME  ( lessRNA_ch )
}

// gwips_RNASeq {
//
//}