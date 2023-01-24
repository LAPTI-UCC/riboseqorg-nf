
process GET_RUN_INFO {
    publishDir "$projectDir/$params.data_dir/runInfos", mode: 'copy'

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

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


"""
TODO: I need to find a way to raise an error if the output file is empty

Sometimes a brokenPipe error leads to the file being empty. 
This appears to be a bug in the ffq package. I had a similar issue with piping to cat in the command line. 
I cannot tell what determines what causes ffq to fail 

"""

process RUN_FFQ {

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        tuple val(GSE),val(srp)

    output:
        path "*.fastq.gz"

    script:
    """
     ffq --ftp $GSE  | jq -r '.[] | .url' | xargs curl -O
    """
}


process WGET_FASTQ {

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        path fastq_url

    output:
        path "*.fastq.gz"

    shell:
    """
    wget -i $fastq_url
    """
}


process FIND_ADAPTERS {
    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*_adpater_report.tsv'

    input:
        file raw_fastq

    output:
        file "${raw_fastq}_adpater_report.tsv"

    script:
        """
        python3 $projectDir/scripts/get_adapters.py -q $raw_fastq -o "${raw_fastq}_adpater_report.tsv"
        """
}


process CLIP_FASTQ {
        
    input:
    file raw_fastq 
    file adapter_report

    output:
    file '*_clipped.fastq' /// into clipped_fastq_channel  ///
	
	script: 
	"""
    cutadapt -j 5 --minimum-length=25 -a "file:$adapter_report" -o $raw_fastq"_clipped.fastq" $raw_fastq
    """
}


process COLLAPSE_FASTQ {
    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*.fastq.gz'

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        path fastq

    output:
        path "*.fastq.gz"

    script:
    """
    python3 $projectDir/scripts/collapse_fastq.py -i $fastq -o ${fastq.baseName}.fastq.gz

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
