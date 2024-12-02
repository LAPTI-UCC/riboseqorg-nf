

process CUTADAPT {
	publishDir "${params.outdir}/trimmed", mode: 'copy'

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

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
