// workflow for processing files for gwips
include { BOWTIE_GENOME } from '../../modules/local/bowtie.nf'
include { SAMTOOLS_INDEX } from '../../modules/local/samtools.nf'
include { GENOME_BAM_TO_BED } from '../../modules/local/bam_to_bed.nf'
include { BEDGRAPH_TO_BIGWIG } from '../../modules/local/bedGraphToBigWig.nf'
include { BEDTOOLS_GENOMECOV } from '../../modules/local/bedtools.nf'

workflow gwips_RiboSeq {

    take:
        lessRNA_ch    

    main:
        genome_sorted_bam_ch    =   BOWTIE_GENOME           ( lessRNA_ch )
	    indexed_bam_ch          =   SAMTOOLS_INDEX   	    ( genome_sorted_bam_ch.genome_sorted_bam )
        bed_ch                  =   BEDTOOLS_GENOMECOV      ( indexed_bam_ch.sorted_bam)
	    // bed_ch                  =   GENOME_BAM_TO_BED       ( indexed_bam_ch.sorted_bam, indexed_bam_ch.sorted_bam_bai, 0 )
		bigwig_ch               =   BEDGRAPH_TO_BIGWIG      ( bed_ch.sorted_beds )
}
