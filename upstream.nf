/* Creating two parameters for the first two arguments of the get_runInfo.py script.
*/

params.ribosome_prof_superset = "/data/ribosome_profiling_superset.csv"
params.data_dir = "data"
project_dir = projectDir 

log.info """\
    R I B O - S E Q    N F    P I P E L I N E
    =========================================
    



"""
.stripIndent()

/// Processes necessary for metadata_flow
include { GET_GSE_REPORT; GET_CSV_FROM_XML; ASSESS_LIBRARY_STRATEGY; CHECK_METADATA_REPORT } from './modules/metadata-tasks.nf'
/// Processes necessary for upstream_flow
include { GET_RUN_INFO; GET_INDIVIDUAL_RUNS; RUN_FFQ; WGET_FASTQ; FIND_ADAPTERS; WRITE_PARAMTERS_YAML } from './modules/upstream-tasks.nf'


/// Workflow to get metadata. Returns a .csv for each GSE with info on their runs.
workflow metadata_flow {

    take: GSE_inputs

    main:
        gse_report_ch = GET_GSE_REPORT( GSE_inputs )
        csv_ch = GET_CSV_FROM_XML( gse_report_ch )
        lib_strat_csv = ASSESS_LIBRARY_STRATEGY ( csv_ch )
        CHECK_METADATA_REPORT ( lib_strat_csv )

}



/// workflow to get the info for the .yaml file.
workflow upstream_flow {

    take: GSE_inputs

    main:
        run_info_ch = GET_RUN_INFO                  ( GSE_inputs )
        runs_ch = GET_INDIVIDUAL_RUNS               ( run_info_ch ) 
        fq_urls_ch = RUN_FFQ                        ( runs_ch.splitText( by: 1 ) ) 

        fastq_path_ch = WGET_FASTQ                  ( fq_urls_ch, GSE_inputs ) 

        // adapter_report_ch = FIND_ADAPTERS           ( WGET_FASTQ.out, GSE_inputs )

    emit:
        run_info_ch
}

workflow yaml_flow {

    take:
        upstream 
        GSE_inputs

    main:
        WRITE_PARAMTERS_YAML    ( upstream, GSE_inputs )

}

workflow {

    ///NB. We could parse through the superset.csv to get the GSEs, instead of relying on this
    GSE_inputs = Channel
        .fromPath("data/ribosome_profiling_superset.csv")
        .splitCsv(header: true)
        .map { row -> tuple("${row.Accession}", "${row.SRA}" )} // use for superset  

        // .map { row -> tuple("${row.GSE}", "${row.SRP}" )} // use for CSV from trips

    main:
        metadata_flow(GSE_inputs)
        // upstream_flow(GSE_inputs)
        // yaml_flow(upstream_flow.out, GSE_inputs)


}
///  IMPORTANT! we need to figure out how to execute processes in the right order: if the yamls is created too son, some info will be missing.
///  We can use .collect() to make sure all the outputs of a process are collected before being set to the following one in the pipeline.

workflow.onComplete{
    println "Pipeline completed at: $workflow.complete"
}
