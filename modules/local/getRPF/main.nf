process GET_RPF {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/getRPF", mode: 'copy'
    publishDir "${params.outdir}/failed_rpf_checks", mode: 'copy', pattern: '*_rpf_checks.txt', saveAs: { filename -> 
        file(filename).text.contains("[FAIL]") ? filename : null 
    }

    conda "bioconda::getrf=1.3.0"  // Assuming getRPF is available in Bioconda

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
    getRPF check ${input_file} \\
        --format collapsed \\
        --output ${prefix}_report.txt \\
        --count-pattern "${count_pattern}" \\
        $args

    # Assuming getRPF writes the RPF check results to a file named ${prefix}_rpf_checks.txt
    # If not, you may need to modify this part based on how getRPF actually outputs this information

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