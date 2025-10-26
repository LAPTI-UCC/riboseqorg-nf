

process RIBOMETRIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${projectDir}/conda/ribometric.yml"

	publishDir "$params.outdir/RiboMetric", mode: 'copy'

    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bam_index)
    path ribometric_annotation

    output:
    tuple val(meta), path("*RiboMetric.html"), emit: html
    tuple val(meta), path("*RiboMetric.json"), emit: json
    tuple val(meta), path("*RiboMetric.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    RiboMetric run \\
        --bam ${transcriptome_bam} \\
        --annotation ${ribometric_annotation} \\
        --threads $task.cpus \\
        --html \\
        --json \\
        --csv \\
        --offset-calculation-method tripsviz \\
        -S 10000000
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_RiboMetric.html
    touch ${prefix}_RiboMetric.json
    touch ${prefix}_RiboMetric.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: 1.0.0
    END_VERSIONS
    """
}