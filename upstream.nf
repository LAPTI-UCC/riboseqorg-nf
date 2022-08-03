/* Creating two parameters for the first two arguments of the get_runInfo.py script.
*/

params.ribosome_prof_superset = "/data/ribosome_profiling_superset.csv"
params.data_dir = "data"
project_dir = projectDir 

include { GET_RUN_INFO; GET_INDIVIDUAL_RUNS; RUN_FFQ; WGET_FASTQ; FIND_ADAPTERS; WRITE_PARAMTERS_YAML } from './modules/upstream-tasks.nf'


workflow {
    gses = ["GSE152556", "GSE112305"]
    GSE_inputs = Channel.fromList(gses)  /* a GSE I want to test. Another candidate is GSE152556*/
    GET_RUN_INFO(GSE_inputs)

    GET_INDIVIDUAL_RUNS(GET_RUN_INFO.out) 
    RUN_FFQ(GET_INDIVIDUAL_RUNS.out) // This will not be the optimal method. I resorted to python because I could not manage I/O with nf or shell 
    WGET_FASTQ(RUN_FFQ.out.flatten()) // This will not be optimal similar to above

    FIND_ADAPTERS(WGET_FASTQ.out)
    WRITE_PARAMTERS_YAML(GET_RUN_INFO.out, FIND_ADAPTERS.out)
}


/*
The way this is structured, we want (and It should) process all elements (the GSEs) given as input.
Such inputs can be emitted from a channel (as a start, try to just give here a list to process)

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

