/* Creating two parameters for the first two arguments of the get_runInfo.py script.
*/

params.ribosome_prof_superset = "/data/ribosome_profiling_superset.csv"
params.data_dir = "data"
project_dir = projectDir 


/// Processes necessary for metadata_flow
include { GET_GSE_REPORT; GET_CSV_FROM_XML; ASSESS_LIBRARY_STRATEGY} from './modules/metadata-tasks.nf'
/// Processes necessary for upstream_flow
include { GET_RUN_INFO; GET_INDIVIDUAL_RUNS; RUN_FFQ; WGET_FASTQ; FIND_ADAPTERS; WRITE_PARAMTERS_YAML } from './modules/upstream-tasks.nf'


/// Workflow to get metadata. Returns a .csv for each GSE with info on their runs.
workflow metadata_flow {

    take: GSE_inputs

    main:
        GET_GSE_REPORT          ( GSE_inputs )
        GET_CSV_FROM_XML        ( GET_GSE_REPORT.out )
        ASSESS_LIBRARY_STRATEGY ( GET_CSV_FROM_XML.out )

}


/// workflow to get the info for the .yaml file.
workflow upstream_flow {

    take: GSE_inputs

    main:
        GET_RUN_INFO            ( GSE_inputs )

        GET_INDIVIDUAL_RUNS     ( GET_RUN_INFO.out ) 
        RUN_FFQ                 ( GET_INDIVIDUAL_RUNS.out ) // This will not be the optimal method. I resorted to python because I could not manage I/O with nf or shell 
        WGET_FASTQ              ( RUN_FFQ.out.flatten(), GSE_inputs ) // This will not be optimal similar to above

        FIND_ADAPTERS           ( WGET_FASTQ.out, GSE_inputs )

    emit:
        GET_RUN_INFO.out
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
        .fromPath("/home/121109636/CSV_reports/GSEs.txt")
        .splitCsv(header: true)
        .map { row -> tuple("${row.Accession}", "${row.SRA}" )}
    
    GSE_inputs.view()
  
    main:
        metadata_flow(GSE_inputs)
        upstream_flow(GSE_inputs)
        yaml_flow(upstream_flow.out, GSE_inputs)


}
///  IMPORTANT! we need to figure out how to execute processes in the right order: if the yamls is created too son, some info will be missing.
///  We can use .collect() to make sure all the outputs of a process are collected before being set to the following one in the pipeline.
