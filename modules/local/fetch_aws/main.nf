process FETCH_RUN {
    tag 'high'

    conda "conda-forge::awscli"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/awscli:1.32.6' :
        'biocontainers/awscli:1.32.6' }"

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        tuple val(study_accession), val(run), val(scientific_name), val(library_type)

    output:
        file "${run}"

    script:
        """
        aws --no-sign-request s3 sync s3://sra-pub-run-odp/sra/${run} .
        """
}