
process GET_RUN_INFO {

    // errorStrategy 'ignore'


    input:
        tuple val(GSE),val(srp)

    output:
        file "*_sraRunInfo.csv"

    script:

        def z = ["4", "5", "6", "7", "8", "9"]
        Random rnd = new Random()

        Sleep_time = (z[rnd.nextInt(z.size)])
        """
        sleep ${Sleep_time}
        esearch -db sra -query ${srp} | efetch -format runinfo -mode text | cat > ${GSE}_sraRunInfo.csv
        """
}


process GET_INDIVIDUAL_RUNS {

    errorStrategy 'ignore'


    input:
        path sraRunInfo

    output:
        // file 'run_*'
        stdout

    script:
    """
    cut -f1 -d, ${sraRunInfo} | tail -n+2 | cat > srr.txt 
    cat srr.txt 
    """
}

process RUN_FFQ {

    errorStrategy 'ignore'

    input:
        val SRR

    output:
        stdout

    script:
    trimmed_string = "${SRR[0..-2]}"
    """
     ffq --ftp $trimmed_string | jq -r .[].url | tr -d '\n'
    """
    }

/// We want all the runs for a study to be published in the same study directory, identified by the GSE

process WGET_FASTQ {
    publishDir "$projectDir/$params.data_dir/$GSE/fastq", mode: 'copy', pattern: '*.fastq.gz'

    errorStrategy 'ignore'

    input:
        path ffq_json
        tuple val(GSE),val(srp)


    output:
        file "*.fastq.gz"

    script:
        def z = ["4", "5", "6", "7", "8", "9"]
        Random rnd = new Random()

        Sleep_time = (z[rnd.nextInt(z.size)])
        """
        sleep ${Sleep_time}
        wget ${url} -P ./
        """
    shell:
    """
    #!/usr/bin/python3

import os
import time
import random

wait_time = random.choice([5,6,7,8,9,10,11,12])
time.sleep(wait_time)

with open('${ffq_json}', 'r') as f:
    url = f.readlines()[0].strip('\\n')
    os.system(f"wget {url} -P ./")

    """

}




process FIND_ADAPTERS {
    publishDir "$projectDir/$params.data_dir/$GSE/fastq", mode: 'copy', pattern: '*_adpater_report.tsv'

    // errorStrategy 'ignore'


    input:
        file raw_fastq
        tuple val(GSE),val(srp)

    output:
        file "${raw_fastq}_adpater_report.tsv"

    script:
    
        """
        python3 $projectDir/scripts/get_adapters.py -q $raw_fastq -o "${raw_fastq}_adpater_report.tsv"
        """
}



process WRITE_PARAMTERS_YAML {

    errorStrategy 'ignore'


    input:
        file sraRunInfo
        tuple val(GSE),val(srp)

    output:
        file "parameters.yaml"

    script:
        """
        python3 $projectDir/scripts/write_parameters_yaml.py -a "$projectDir/$params.data_dir/$GSE/fastq" -s $projectDir/annotation_inventory/annotation_inventory.sqlite -r $sraRunInfo -o parameters.yaml
        """
}
