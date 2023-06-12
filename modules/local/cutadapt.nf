

process CUTADAPT {
	publishDir "${params.study_dir}/trimmed", mode: 'copy'

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
