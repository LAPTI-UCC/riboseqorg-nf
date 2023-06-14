process FETCH_RUN {
    tag 'high'

    publishDir "$projectDir/$params.study_dir/fastq", mode: 'copy', pattern: '*.fastq.gz'

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