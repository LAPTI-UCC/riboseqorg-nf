import java.util.Random

process FASTQ_DL {
    tag 'medium'

    conda '${projectDir}/conda/fastq-dl.yaml'

    publishDir "$projectDir/$params.study_dir/fastq", mode: 'copy', pattern: '*.fastq.gz'

    input:
        tuple val(study_accession), val(run), val(scientific_name), val(library_type)

    output:
        file "*.fastq.gz"

    script:
        def sleepDuration = random.nextInt(10) + 1

        sleep(sleepDuration)

        """
        fastq-dl -a $run --cpus 5
        """
}