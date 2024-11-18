

process FASTQC {

    tag 'medium'

    conda "${projectDir}/conda/fastqc.yaml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

	publishDir "${params.study_dir}/fastqc", mode: 'copy'
	
	input:
	    file fastq 

	output:
	    path "*_fastqc.html", emit: fastqc_html
        path "*/fastqc_data.txt", emit: fastqc_data

    script:
        """
        fastqc --extract -q $fastq --adapters $projectDir/scripts/adapter_list.tsv
        """
}
