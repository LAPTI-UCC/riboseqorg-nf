

process FASTERQ_DUMP {
    tag 'high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:2f4a4c900edd6801ff0068c2b3048b4459d119eb-0' :
        'biocontainers/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:2f4a4c900edd6801ff0068c2b3048b4459d119eb-0' }"

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        file run

    output:
        file "*.fastq"

    script:
        """
        fasterq-dump -f -p -e $task.cpus --split-files --skip-technical ${run}
        """
}