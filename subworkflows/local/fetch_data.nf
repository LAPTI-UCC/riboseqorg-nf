// workflow for fetching data

include { LOCATE } from '../../modules/local/locate/main.nf'
include { collapse } from '../../subworkflows/local/collapse.nf'

workflow fetch_data {

    take: samples_ch

    main:
        if (params.fetch == null) {
            error "params.fetch must be defined as true or false"
        }
        
        LOCATE( samples_ch)

        // Handle files that need processing
        LOCATE.out.needs_processing
            .set { needs_processing }
        
        // Run the collapse workflow for runs that need processing, only if fetch is true
        if (params.fetch) {
            collapse(needs_processing)
            newly_collapsed_reads = collapse.out.collapsed_fastq
        } else {
            newly_collapsed_reads = Channel.empty()
            needs_processing.view { run -> 
                log.warn "Run $run needs processing but fetch is set to false. This run will be skipped."
            }
        }

        // Combine existing collapsed reads with newly processed ones
        collapsed_reads = LOCATE.out.collapsed_reads
            .mix(newly_collapsed_reads)
            .groupTuple()
            .map { run, files -> tuple(run, files[0]) }  // Take the first file if there are duplicates
            
    emit:
        collapsed_reads
}
