

process CUTADAPT {
	publishDir "${params.outdir}/trimmed", mode: 'copy'

    conda "bioconda::cutadapt=5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:5.2--py39hbcbf7aa_0' :
        'biocontainers/cutadapt:5.2--py39hbcbf7aa_0' }"

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
