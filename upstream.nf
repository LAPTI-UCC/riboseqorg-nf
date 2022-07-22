/* Creating two parameters for the first two arguments of the get_runInfo.py script.
*/

params.ribosome_prof_superset = "/data/ribosome_profiling_superset.csv"
params.data_folder = "/data"
project_dir = projectDir 

params.study_dir = "/home/115316376"




process GET_RUN_INFO {

    input:
        val GSE

    output:
        file "${GSE}_sraRunInfo.csv"

    script:
        """
        python3 $project_dir/scripts/get_runInfo.py $project_dir${params.ribosome_prof_superset} $project_dir${params.data_folder} $GSE ${GSE}_sraRunInfo.csv
        """
}

process GET_INDIVIDUAL_RUN_INFOS {
    publishDir "/home/115316376/individual_runInfos", mode: 'copy', pattern: '*_sraRunInfo.csv'

    input:
        file sraRunInfo

    output:
        file '*_sraRunInfo.csv'

    script:
        """
        python3 $project_dir/scripts/split_runInfo_to_rows.py -r $sraRunInfo -o ./
        """
}

process GET_FASTQ {

    input:
        path sraRunInfo

    output:
        file '*.fastq.gz'

    script:
        """
        python3 $project_dir/scripts/ffq_fetch_fastq.py -r $sraRunInfo -o ./
        """
}


process FIND_ADAPTERS {
    publishDir "$params.study_dir/adapter_reports", mode: 'copy', pattern: '*_adpater_report.fa'


    input:
        file raw_fastq

    output:
        file "${raw_fastq}_adapter_report.fa"

    script:
        """
        python3 $project_dir/scripts/get_adapters.py -q $raw_fastq -o "${raw_fastq}_adapter_report.fa"
        """
}



process WRITE_PARAMTERS_YAML {


    input:
        file sraRunInfo
        path find_adapters

    output:
        file "paramters.yaml"

    script:
        """
        python3 $project_dir/scripts/write_paramaters_yaml.py -a "${params.study_dir}/adapter_reports" -s $project_dir/annotation_inventory/annotation_inventory.sqlite -r $sraRunInfo -o paramaters.yaml
        """
}


workflow {
    GSE_inputs = Channel.of("GSE112305")  /* a GSE I want to test. Another candidate is GSE152556*/
    GET_RUN_INFO(GSE_inputs)
    GET_INDIVIDUAL_RUN_INFOS(GET_RUN_INFO.out) /* this outputs a string of filenames and I want a channel */
    GET_FASTQ(GET_INDIVIDUAL_RUN_INFOS.out.flatten())
    FIND_ADAPTERS(GET_FASTQ.out)
    // WRITE_PARAMTERS_YAML(GET_RUN_INFO.out, FIND_ADAPTERS.out)
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

