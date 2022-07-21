/* Creating two parameters for the first two arguments of the get_runInfo.py script.
*/

params.ribosome_prof_superset = "/data/ribosome_profiling_superset.csv"
params.data_folder = "/data"
project_dir = projectDir 

process GET_RUN_INFO {

input:
val GSE

output:
file "${GSE}_sraRunInfo.csv"

script:
"""
python3 $project_dir/scripts/get_runInfo.py $project_dir${params.ribosome_prof_superset} $project_dir${params.data_folder} $GSE $GSE"_sraRunInfo.csv"
"""
}


workflow {
GSE_inputs = Channel.of("GSE158141")  /* a GSE I want to test. Another candidate is GSE152556*/
GET_RUN_INFO(GSE_inputs)
GET_RUN_INFO.out.view()
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

