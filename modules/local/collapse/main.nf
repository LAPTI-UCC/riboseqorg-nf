process COLLAPSE_FASTQ {
    tag "$meta.id"

    conda "${projectDir}/conda/RDP-tools.yml"

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rdp-tools:1.0.0--pyhdfd78af_0' :
        'biocontainers/rdp-tools:1.0.0--pyhdfd78af_0' }"

    publishDir "${params.outdir}/collapsed_fa", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*collapsed.fa.gz"), emit: collapsed_fastq

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    RDP-Tools collapse $args $fastq
    gzip *fastq.collapsed.fa
    """ 

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_collapsed.fa.gz
    """
}