

include { STAR_ALIGN } from '../../modules/local/STAR/main'
include { SAMTOOLS_COORD_SORT } from '../../modules/local/samtools/samtools_coord_sort/main'
include { SAMTOOLS_INDEX } from '../../modules/local/samtools/samtools_index/main'

workflow ALIGNMENT {
    take:
    samples
    star_index
    gtf

    main:
    STAR_ALIGN(samples, star_index, gtf)
    SAMTOOLS_COORD_SORT(STAR_ALIGN.out.transcriptome_bam)
    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)

    emit:
    bams = SAMTOOLS_COORD_SORT.out.bam
    bai = SAMTOOLS_INDEX.out.bai
    genome_bam = STAR_ALIGN.out.bam

}