// Just works for fastqc output for now

process MULTIQC {

	publishDir "$params.study_dir/multiqc", mode: 'copy'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.22.3--pyhdfd78af_0' :
        'biocontainers/multiqc:1.22.3--pyhdfd78af_0' }"

	input:
	file ('fastqc/*')

	output:
	file "multiqc_report.html"
	
	"""
	multiqc $params.study_dir/fastqc
	"""
}