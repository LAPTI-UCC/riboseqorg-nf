include { SAMTOOLS_NAME_SORT } from '../../modules/local/samtools/samtools_name_sort/main'
include { BAM_TO_SQLITE } from '../../modules/local/bam_to_sqlite/main'
include { EXTRACT_OFFSETS } from '../../modules/local/extract_offsets/main'
include { BAM_TO_BED } from '../../modules/local/bam_to_bed/main'
include { BEDGRAPH_TO_BIGWIG } from '../../modules/local/bedGraphToBigWig/main'

workflow POST_PROCESSING {
    take:
    bams 
    genome_bam
    bai  
    annotation_sqlite
    chrom_sizes

    main:
    // Combine sorted BAMs with their index files
    sorted_bams_with_index = genome_bam.join(bai)

    // Name sort for BAM_TO_SQLITE (only this branch needs name sorting)
    SAMTOOLS_NAME_SORT(bams)
    BAM_TO_SQLITE(SAMTOOLS_NAME_SORT.out.bam, annotation_sqlite)
    EXTRACT_OFFSETS(BAM_TO_SQLITE.out.sqlite)

    // Join coordinate-sorted BAMs (with index) and offsets
    bam_bai_and_offsets = sorted_bams_with_index
        .join(EXTRACT_OFFSETS.out.offsets, by: [0])

    // Run BAM_TO_BED
    BAM_TO_BED(bam_bai_and_offsets)

    unpacked_bedgraphs = BAM_TO_BED.out.bedgraph
        .flatMap { meta, bedgraphs ->
            bedgraphs.collect { bedgraph ->
                [meta, bedgraph]
            }
        }

    BEDGRAPH_TO_BIGWIG(unpacked_bedgraphs, chrom_sizes)


    emit:
    bedgraphs = BAM_TO_BED.out.bedgraph
    coord_sorted_bams = bams
    bam_indices = bai
    sqlite_dbs = BAM_TO_SQLITE.out.sqlite
    offsets = EXTRACT_OFFSETS.out.offsets
}