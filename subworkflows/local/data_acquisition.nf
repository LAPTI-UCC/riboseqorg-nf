

include { LOCATE } from '../../modules/local/locate/main.nf'
include { collapse } from './collapse.nf'

workflow DATA_ACQUISITION {
    take:
    sample_sheet

    main:
    // Validate required parameters
    if (params.fetch == null) {
        error "params.fetch must be defined as true or false"
    }
    
    // Parse sample sheet and create samples channel
    samples_ch = Channel
        .fromPath(sample_sheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            def meta = [
                id: row.Run,
                study_accession: row.study_accession,
            ]
            [ meta, row.Run ]
        }
    
    // Locate existing collapsed reads or runs that need processing
    LOCATE(samples_ch)
    LOCATE.out.needs_processing
        .map { meta, id, file -> id }
        .collectFile(name: "${params.outdir}/runs_needing_processing.txt", newLine: true)
        
    // Handle different combinations of fetch and force_fetch
    if (params.fetch) {
        if (params.force_fetch) {
            // fetch = true, force_fetch = true
            // Process all samples, including existing collapsed reads
            needs_processing = LOCATE.out.needs_processing
                .mix(
                    LOCATE.out.collapsed_reads.map { meta, collapsed_file -> 
                        // Use the existing file path for force fetch
                        [meta, meta.id, collapsed_file]
                    }
                )
            collapsed_reads = Channel.empty()
        } else {
            // fetch = true, force_fetch = false
            // Process only samples that need processing
            needs_processing = LOCATE.out.needs_processing
            collapsed_reads = LOCATE.out.collapsed_reads
        }
        
        // Run the collapse workflow
        collapse(needs_processing)
        newly_collapsed_reads = collapse.out.collapsed_fastq
    } else {
        if (params.force_fetch) {
            // fetch = false, force_fetch = true
            // This is a contradictory state, we should warn the user
            log.warn "fetch is set to false but force_fetch is true. This combination is not valid. No data will be fetched or processed."
        } else {
            // fetch = false, force_fetch = false
            // Don't fetch or process any data
            log.info "fetch is set to false. No data will be fetched or processed."
        }
        
        needs_processing = Channel.empty()
        newly_collapsed_reads = Channel.empty()
        collapsed_reads = LOCATE.out.collapsed_reads
    }

    // Combine existing collapsed reads with newly processed ones
    all_collapsed_reads = collapsed_reads
        .mix(newly_collapsed_reads)
        .groupTuple()
        .map { run, files -> tuple(run, files[0]) }  // Take the first file if there are duplicates
        
    emit:
    samples = all_collapsed_reads
}