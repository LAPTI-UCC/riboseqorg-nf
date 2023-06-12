// workflow for processing files for trips
include { BOWTIE_TRANSCRIPTOME } from '../../modules/local/bowtie.nf'


workflow trips_RiboSeq {


    take: 
        lessRNA_ch    

    main:
        transcriptome_bam_ch    =   BOWTIE_TRANSCRIPTOME  ( lessRNA_ch )

    emit:
        transcriptome_bam_ch.transcriptome_bam
}

// trips_RNASeq {
//
//}