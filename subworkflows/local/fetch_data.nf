// workflow for fetching data

/// Processes necessary for upstream_flow
include { GET_RUN_INFO; CLIP_FASTQ; COLLAPSE_FASTQ; FIND_ADAPTERS; WRITE_PARAMTERS_YAML; FASTQ_DL            } from './modules/upstream-tasks.nf'

include { FASTQC } from './modules/processing-tasks.nf'


workflow fetch_data {

    take: GSE_inputs

    main:
        run_info_ch         = GET_RUN_INFO            ( GSE_inputs )
        runs_ch             = run_info_ch.splitCsv(header: true).map { row -> tuple("${row.Run}", params.GSE )}
        fastq_path_ch       = FASTQ_DL                ( runs_ch )

    emit:
        fastq_path_ch
}
