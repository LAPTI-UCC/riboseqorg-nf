
process GET_RUN_INFO {

    // errorStrategy 'ignore'
    publishDir "$projectDir/$params.data_dir/runInfos", mode: 'copy'

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




process RUN_FFQ {

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        tuple val(GSE),val(srp)

    output:
        // file "*.json"
        path "*.txt"

    script:
    """
     ffq --ftp $GSE | cat > ${GSE}.json
     cat ${GSE}.json | jq -r .[].url > outfile.txt

    """
}

process GET_URL {

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        path json
        // tuple val(GSE),val(srp)

    output:
        path "*.txt"

    shell:
    """
    #!/bin/python3

    import json

    with open('${json}', 'r', encoding='utf-8') as jsonfile:
        jsonfile.seek(0)
        data = json.load(jsonfile)

        url = data[0]["url"]

    outfile = open("outfile.txt", 'w')

    outfile.write(url)
    outfile.close()
    """

}


process WGET_FASTQ {
    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*.fastq.gz'

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        path fastq_url
        // tuple val(GSE),val(srp)


    output:
        path "*.fastq.gz"

    script:
        """
        wget -i $fastq_url
        """

}




process FIND_ADAPTERS {
    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*_adpater_report.tsv'

    // errorStrategy 'ignore'


    input:
        file raw_fastq
        // tuple val(GSE),val(srp)

    output:
        file "${raw_fastq}_adpater_report.tsv"

    script:
    
        """
        python3 $projectDir/scripts/get_adapters.py -q $raw_fastq -o "${raw_fastq}_adpater_report.tsv"
        """
}



process WRITE_PARAMTERS_YAML {
    publishDir "$projectDir/$params.data_dir/$params.GSE", mode: 'copy', pattern: '*.yaml'

    errorStrategy 'ignore'

    input:
        path adapter_report_ch
        file sraRunInfo
        tuple val(GSE),val(srp)

    output:
        file "parameters.yaml"

    script:
        """
        python3 $projectDir/scripts/write_parameters_yaml.py -a "$projectDir/$params.data_dir/$GSE/fastq" -s $projectDir/annotation_inventory/annotation_inventory.sqlite -r $sraRunInfo -o parameters.yaml
        """
}
