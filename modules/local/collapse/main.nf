process COLLAPSE_FASTQ {
    tag "$meta.id"

    conda "${projectDir}/conda/RDP-tools.yml"
    container "ghcr.io/lapti-ucc/riboseqorg-nf-rdp-tools:latest"

    label 'process_high'

    publishDir "${params.outdir}/collapsed_fa", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*collapsed.fa"), emit: collapsed_fasta

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    RDP-Tools collapse $args $fastq
    """ 

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_collapsed.fa
    """
}