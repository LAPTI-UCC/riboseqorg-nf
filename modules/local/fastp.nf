process FASTP {
    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*.fastq.gz'

    input:
    file raw_fastq 
    file adapter_report

    output:
    file '*_clipped.fastq' 
	
	script: 
	"""
    fastp --in $raw_fastq --out $raw_fastq.baseName + '_clipped.fastq' --length_required 15 --adapter_sequence file:$adapter_report
    """
}
