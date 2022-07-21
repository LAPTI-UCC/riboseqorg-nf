/*
params.smthing = "path_to_csv_superset.csv"  --> this ideally should be an argument given to nextflow at the start from terminal
params.smthingelse = "path?to/data/folder"   --> as above

process GET_RUN_INFO {

    input:
    val GSE
    output:
    path CSV, emit: csv
    path H_CSV, emit: header_csv     ----> This relates to the problem to define the input for process 2. 
    
""""
python3 ./scripts/get_runInfo.py $params.smthing $params.smthingelse $GSE
""""
}
*/
/* 
workflow {
GSE_inputs = Channel.from( GSE , GSE , GSE , GSE )
GET_RUN_INFO(GSE_inputs)
}
*/
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

