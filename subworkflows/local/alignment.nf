

include { STAR_ALIGN            } from '../../modules/local/STAR/main'
include { BOWTIE_ALIGN_SORT     } from '../../modules/local/bowtie/align/main'
include { BOWTIE_REMOVE         } from '../../modules/local/bowtie/remove/main'
include { RIBODETECTOR_REMOVE   } from '../../modules/local/ribodetector/remove/main'
include { SAMTOOLS_COORD_SORT   } from '../../modules/local/samtools/samtools_coord_sort/main'
include { SAMTOOLS_INDEX        } from '../../modules/local/samtools/samtools_index/main'

workflow ALIGNMENT {
    take:
    samples
    star_index
    bowtie_index
    rRNA_index
    gtf

    main:
    STAR_ALIGN          (samples, star_index, gtf)

    // Choose rRNA removal method
    Channel
        .of( params.rrna_removal_method ?: 'bowtie' )
        .first()
        .view{ "Using rRNA removal method: ${it}" }

    if ((params.rrna_removal_method ?: 'bowtie') == 'ribodetector') {
        RIBODETECTOR_REMOVE (samples, rRNA_index)
        BOWTIE_ALIGN_SORT   (RIBODETECTOR_REMOVE.out.fastq, bowtie_index)
    } else {
        BOWTIE_REMOVE       (samples, rRNA_index)
        BOWTIE_ALIGN_SORT   (BOWTIE_REMOVE.out.fastq, bowtie_index)
    }

    emit:
    transcriptome_bam   = BOWTIE_ALIGN_SORT.out.bam
    genome_bam          = STAR_ALIGN.out.bam
}
