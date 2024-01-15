process FASTERQ_DUMP {
    tag 'high'

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