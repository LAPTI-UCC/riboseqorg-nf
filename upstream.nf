/* Creating two parameters for the first two arguments of the get_runInfo.py script.
*/
nextflow.enable.dsl=2


params.ribosome_prof_superset = "/data/ribosome_profiling_superset.csv"
params.data_dir = "data"
project_dir = projectDir 

log.info """\
    R I B O - S E Q    N F    P I P E L I N E
    =========================================
    

"""


/// Processes necessary for upstream_flow
include { GET_RUN_INFO; RUN_FFQ; CLIP_FASTQ; COLLAPSE_FASTQ; WGET_FASTQ; FIND_ADAPTERS; WRITE_PARAMTERS_YAML            } from './modules/upstream-tasks.nf'



workflow upstream_flow {

    take: GSE_inputs

    main:
        run_info_ch         = GET_RUN_INFO            ( GSE_inputs )
        runs_ch             = run_info_ch.splitCsv(header: true).map { row -> tuple("${row.Run}", params.GSE )}
        fastq_path_ch       = RUN_FFQ                 ( runs_ch ) 
        adapter_report_ch   = FIND_ADAPTERS           ( fastq_path_ch )
        clipped_ch          = CLIP_FASTQ              ( fastq_path_ch, adapter_report_ch )
        collapsed_ch        = COLLAPSE_FASTQ          ( clipped_ch )
        params_ch           = WRITE_PARAMTERS_YAML    (adapter_report_ch.collect(), run_info_ch, GSE_inputs ) 

    emit:
        params_ch
}


workflow {

    GSE_inputs = Channel
        .fromPath("data/GSE69921_test.csv")
        .splitCsv(header: true)
        .map { row -> tuple("${row.Accession}", "${row.SRA}" )} // use for superset  

        // .map { row -> tuple("${row.GSE}", "${row.SRP}" )} // use for CSV from trips

    main:
        upstream_flow( GSE_inputs )


}

workflow.onComplete{
    println "Pipeline completed at: $workflow.complete"
}
