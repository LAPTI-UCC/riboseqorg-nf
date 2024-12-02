process MULTIQC {
    conda "conda-forge::python=3.11 conda-forge::typing_extensions bioconda::multiqc=1.17.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.22.3--pyhdfd78af_0' :
        'biocontainers/multiqc:1.22.3--pyhdfd78af_0' }"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path "*"

    output:
    path("*multiqc_report.html"), emit: report
    path("*_data"), emit: data
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    multiqc \\
        --force \\
        $args \\
        -o . \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:    """
    touch stub_multiqc_report.html
    mkdir stub_data
    touch stub_data/multiqc_general_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.22.3
    END_VERSIONS
    """
}