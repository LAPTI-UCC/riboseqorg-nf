
process GET_RUN_INFO {

    input:
        val GSE

    output:
        path "${GSE}_sraRunInfo.csv"

    script:
        """
        python3 $projectDir/scripts/get_runInfo.py $projectDir/${params.ribosome_prof_superset} $projectDir/${params.data_dir} $GSE ${GSE}_sraRunInfo.csv
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
        print(line)
        line = line.strip('\\n')
        os.system(f"ffq --ftp {line} | jq -r .[].url | cat > ./{line}.json")


    """
    }


process WGET_FASTQ {
    publishDir "../$params.data_dir/adapter_reports", mode: 'copy', pattern: '*_adpater_report.fa'

    input:
        path ffq_json

    output:
        file "*.fastq.gz"

    shell:
    """
    #!/usr/bin/python3

import os
import time
import random

wait_time = random.choice([5,15,25,35,45,55,60])
time.sleep(wait_time)

with open('${ffq_json}', 'r') as f:
    url = f.readlines()[0].strip('\\n')
    os.system(f"wget {url} -P ./")

    """

}




process FIND_ADAPTERS {
    publishDir "../$params.data_dir/adapter_reports", mode: 'copy', pattern: '*_adpater_report.fa'


    input:
        file raw_fastq

    output:
        file "${raw_fastq}.fa"

    script:
    
        """
        python3 ../scripts/get_adapters.py -q $raw_fastq -o "${raw_fastq}.fa"
        """
}



process WRITE_PARAMTERS_YAML {


    input:
        file sraRunInfo
        path find_adapters

    output:
        file "parameters.yaml"

    script:
        """
        python3 ../scripts/write_parameters_yaml.py -a "$project_dir/$params.data_dir/$find_adapters.simpleName/adapter_reports" -s $project_dir/annotation_inventory/annotation_inventory.sqlite -r $sraRunInfo -o parameters.yaml
        """
}
