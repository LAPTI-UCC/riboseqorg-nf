process GWIPS_INSERTS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9.5 conda-forge::sqlite=3.39.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9.5' :
        'quay.io/biocontainers/python:3.9.5' }"

    input:
    tuple val(meta), path(run_metadata)
    path study_metadata
    path annotation_inventory_sqlite

    output:
    tuple val(meta), path("*.txt"), emit: inserts
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 $projectDir/scripts/write_GWIPS_inserts.py \\
        -s ${study_metadata} \\
        -m ${run_metadata} \\
        --db ${annotation_inventory_sqlite} \\
        --output ${prefix}_gwips_inserts.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sqlite: \$(sqlite3 --version | sed 's/^.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_gwips_inserts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.5
        sqlite: 3.39.3
    END_VERSIONS
    """
}