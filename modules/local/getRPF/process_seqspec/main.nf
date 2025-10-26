process PROCESS_SEQSPEC {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pyyaml=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    publishDir "${params.outdir}/seqspec_adapters", mode: 'copy'

    input:
    tuple val(meta), path(seqspec_file)

    output:
    tuple val(meta), path("*.adapters.fa"), emit: adapter_fasta
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    process_seqspec.py \\
        --seqspec ${seqspec_file} \\
        --output-fasta ${prefix}.adapters.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.adapters.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
