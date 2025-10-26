

process EXTRACT_OFFSETS {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9 conda-forge::sqlite=3.39.3"

    publishDir "$params.outdir/offsets", mode: 'copy'

    input:
    tuple val(meta), path(sqlite_file)

    output:
    tuple val(meta), path("*offsets.txt"), emit: offsets
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python $projectDir/bin/extract_offsets.py \\
        -d ${sqlite_file} \\
        -o ${prefix}.offsets.txt \\
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
    touch ${prefix}.offsets.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.5
        sqlite: 3.39.3
    END_VERSIONS
    """
}