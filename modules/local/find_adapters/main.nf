process FIND_ADAPTERS {
    tag "$meta.id"

    conda "conda-forge::python=3.9 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    publishDir "${params.outdir}/adapter_reports", mode: 'copy'
    
    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
    tuple val(meta), path(raw_fastq)
    tuple val(meta2), path(fastqc_data)

    output:
    tuple val(meta), path("*_adapter_report.fa"), emit: adapter_report
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 $projectDir/bin/get_adapters.py \\
        -i $fastqc_data \\
        -a $projectDir/resources/adapter_list.tsv \\
        -o "${prefix}_adapter_report.fa" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_adapters.py: \$(python3 $projectDir/bin/get_adapters.py --version 2>&1 | sed 's/get_adapters.py v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_adapter_report.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
        get_adapters.py: 1.0.0
    END_VERSIONS
    """
}