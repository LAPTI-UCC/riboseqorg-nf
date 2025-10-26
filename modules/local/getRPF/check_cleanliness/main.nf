process CHECK_CLEANLINESS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${projectDir}/conda/getRPF.yml"

    publishDir "${params.outdir}/getRPF/check", mode: 'copy'

    // errorStrategy 'ignore' 

    input:
    tuple val(meta), path(input_file)
    val(count_pattern)

    output:
    tuple val(meta), path("*_report.txt"), emit: report
    tuple val(meta), path("*rpf_checks.txt"), emit: rpf_checks
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    getRPF check-cleanliness ${input_file} \\
        --format collapsed \\
        --output . \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getrf: \$(getRPF --version 2>&1 | sed 's/getRPF version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.txt
    touch ${prefix}_rpf_checks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getrf: 1.3.0
    END_VERSIONS
    """
}