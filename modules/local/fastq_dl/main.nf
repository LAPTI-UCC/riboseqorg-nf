process FASTQ_DL {
    tag "${meta.id}"

    conda "bioconda::fastq-dl=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-dl:2.0.1--pyhdfd78af_0' :
        'biocontainers/fastq-dl:2.0.1--pyhdfd78af_0' }"

    publishDir "${params.outdir}/fastq", mode: 'copy', pattern: '*.fastq.gz'

    errorStrategy 'ignore'

    input:
    tuple val(meta), val(run), path(needs_processing)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastq-dl \\
        -a $run \\
        --cpus $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-dl: \$(fastq-dl --version | sed 's/fastq-dl //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_1.fastq.gz
    touch ${prefix}_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-dl: 2.0.1
    END_VERSIONS
    """
}