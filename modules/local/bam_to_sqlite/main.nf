

process BAM_TO_SQLITE {
    tag "${meta.id}"
    label 'process_medium'

	publishDir "$params.outdir/sqlites", mode: 'copy'

    conda "bioconda::pysam=0.23.3 conda-forge::sqlite=3.39.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369acef7665b99911ea8ddb73d4-0' :
        'quay.io/biocontainers/mulled-v2-9adca5a7d3b24119897cfc20386da6c7fa47bdab:77c1885b47edc369acef7665b99911ea8ddb73d4-0' }"

    input:
    tuple val(meta), path(sorted_bam)
    path annotation_sqlite

    output:
    tuple val(meta), path("*.sqlite"), emit: sqlite
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 $projectDir/bin/bam_to_sqlite.py \\
        --bam ${sorted_bam} \\
        --annotation ${annotation_sqlite} \\
        --output ${prefix}.sqlite \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        sqlite: \$(sqlite3 --version | sed 's/^.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sqlite

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.5
        pysam: 0.19.0
        sqlite: 3.39.3
    END_VERSIONS
    """
}