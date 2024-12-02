process FASTP {
    tag "$meta.id"

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.json'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.html'


    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(adapter_fasta)

    output:
    tuple val(meta), path("*_clipped.fastq.gz"), emit: trimmed_fastq
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastp \\
        -i $reads \\
        -o ${prefix}_clipped.fastq.gz \\
        $args \\
        --length_required 20 \\
        --adapter_fasta $adapter_fasta \\
        --json ${prefix}_fastp.json \\
        --html ${prefix}_fastp.html \\
        --thread $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clipped.fastq.gz
    touch ${prefix}.json
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: 0.23.4
    END_VERSIONS
    """
}