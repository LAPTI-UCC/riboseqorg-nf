

process CUTADAPT {
    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*.fastq.gz'

    input:
    file raw_fastq 
    file adapter_report

    output:
    file '*_clipped.fastq' 
	
	script: 
	"""
    cutadapt -j $task.cpus --minimum-length=15 -a "file:$adapter_report" -o $raw_fastq"_clipped.fastq" $raw_fastq
    """
}
