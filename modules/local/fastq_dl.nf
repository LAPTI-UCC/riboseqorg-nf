import java.util.Random

process FASTQ_DL {
    tag 'medium'

    publishDir "$projectDir/$params.study_dir/fastq", mode: 'copy', pattern: '*.fastq.gz'

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        tuple val(study_accession), val(run), val(scientific_name), val(library_type)

    output:
        file "*.fastq.gz"

    script:
        def sleepDuration = random.nextInt(100) + 4

        sleep(sleepDuration)

        """
        fastq-dl -a $run --cpus $task.cpus --silent
        """
}