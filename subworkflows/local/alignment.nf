

include { STAR_ALIGN } from '../../modules/local/STAR/main'
include { BOWTIE_ALIGN_SORT } from '../../modules/local/bowtie/align/main'
include { SAMTOOLS_COORD_SORT } from '../../modules/local/samtools/samtools_coord_sort/main'
include { SAMTOOLS_INDEX } from '../../modules/local/samtools/samtools_index/main'

workflow ALIGNMENT {
    take:
    samples
    star_index
    bowtie_index
    gtf

    main:
    STAR_ALIGN(samples, star_index, gtf)
    BOWTIE_ALIGN_SORT(samples, bowtie_index)
    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)

    emit:
    transcriptome_bam = BOWTIE_ALIGN_SORT.out.bam
    genome_bam = STAR_ALIGN.out.bam

}