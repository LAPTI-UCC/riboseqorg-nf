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
        WGET_FASTQ              ( RUN_FFQ.out.flatten() ) // This will not be optimal similar to above

        FIND_ADAPTERS           ( WGET_FASTQ.out )
        WRITE_PARAMTERS_YAML    ( GET_RUN_INFO.out, FIND_ADAPTERS.out )

}

workflow {

    ///NB. We could parse through the superset.csv to get the GSEs, instead of relying on this
    GSE_inputs = Channel
        .fromPath("/home/121109636/CSV_reports/GSEs.txt")
        .splitCsv(header: true)
        .map { row -> tuple("${row.Accession}", "${row.title}" }
    
    GSE_inputs.view()
  
    main:
        metadata_flow(GSE_inputs)
        upstream_flow(GSE_inputs)


}
/*

Now, there are several aspects to consider:
1) Is the script "get_runInfo.py" still publishing the directories for each study with the csv inside?
    2) If yes, can we name a "with_header.csv" file as output to feed to the next process?
    
    3) If not, we need to find a way with PublishDir (but it seems unlikely. The script should still be executed in all of its parts)

Strategy: write a one-process workflow, and see the answer to these questions. Once you see that, you can reason on how to link the preprocessing script to it.

Ideally, the output from one process should be the input to the next...I mean, we want this to work "vertically":
a GSE is taken to produce a .csv --> the .csv is used to create a yaml.

we don,t want to produce all csv and then work through each.
to that extent, it is vital to have the input of process 2 to be the output of process 1.

I think I can use emit to that extent, but still I am not sure of the final result.

*/

