// workflow for processing files for trips
include { BOWTIE_TRANSCRIPTOME } from '../../modules/local/bowtie.nf'
include { BAM_TO_SQLITE } from '../../modules/local/riboseqorg.nf'
include { RIBOMETRIC } from '../../modules/local/ribometric.nf'
include { SAMTOOLS_NAME_SORT; SAMTOOLS_INDEX; SAMTOOLS_COORD_SORT } from '../../modules/local/samtools.nf'


workflow trips_RiboSeq {

    take: 
        lessRNA_ch    

    main:
        transcriptome_bam_ch    =   BOWTIE_TRANSCRIPTOME  ( lessRNA_ch )
        sorted_bam_ch           =   SAMTOOLS_COORD_SORT    ( transcriptome_bam_ch.transcriptome_bam )
        indexed_bam_ch          =   SAMTOOLS_INDEX   	  ( sorted_bam_ch )
        ribometric_reports_ch   =   RIBOMETRIC            ( indexed_bam_ch )
        name_sorted_ch          =   SAMTOOLS_NAME_SORT    ( transcriptome_bam_ch.transcriptome_bam )

        sqlite_ch               =   BAM_TO_SQLITE         ( name_sorted_ch )
}

// trips_RNASeq {
//
//}