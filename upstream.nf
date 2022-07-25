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
        path "${GSE}_sraRunInfo.csv"

    script:
        """
        python3 $project_dir/scripts/get_runInfo.py $project_dir${params.ribosome_prof_superset} $project_dir${params.data_folder} $GSE ${GSE}_sraRunInfo.csv
        """
}


process GET_INDIVIDUAL_RUNS {

    input:
        path sraRunInfo

    output:
        file '*.txt'

    script:
    """
    cut -f1 -d, ${sraRunInfo} | tail -n+2 | cat > srrs.txt
    """
}

process RUN_FFQ {

    input:
        file SRR

    output:
        file "*.json"

    shell:
    """
    #!/usr/bin/python3

import os

with open('${SRR}', 'r') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip('\\n')
        os.system(f"ffq --ftp {line} | jq -r .[].url | cat > ./{line}.json")


    """
    }


process WGET_FASTQ {

    input:
        path ffq_json

    output:
        file "*.fastq.gz"

    shell:
    """
    #!/usr/bin/python3

import os

with open('${ffq_json}', 'r') as f:
    url = f.readlines()[0].strip('\\n')
    os.system(f"wget {url} -P ./")

    """

}


process WGET_FASTQ_SHELL {

    input:
        path ffq_json

    output:
        file "*.fastq.gz"

    script:
    """
    cat ${ffq_json} | wget 
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
        python3 $project_dir/scripts/write_parameters_yaml.py -a "${params.study_dir}/adapter_reports" -s $project_dir/annotation_inventory/annotation_inventory.sqlite -r $sraRunInfo -o paramaters.yaml
        """
}


workflow {
    GSE_inputs = Channel.of("GSE112305")  /* a GSE I want to test. Another candidate is GSE152556*/
    GET_RUN_INFO(GSE_inputs)

    GET_INDIVIDUAL_RUNS(GET_RUN_INFO.out) 
    GET_INDIVIDUAL_RUNS.out.view()
    RUN_FFQ(GET_INDIVIDUAL_RUNS.out) // This will not be the optimal method. I resorted to python because I could not manage I/O with nf or shell 
    WGET_FASTQ_SHELL(RUN_FFQ.out.flatten())
    // WGET_FASTQ(RUN_FFQ.out.flatten()) // This will not be optimal similar to above

    // FIND_ADAPTERS(WGET_FASTQ.out)
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

