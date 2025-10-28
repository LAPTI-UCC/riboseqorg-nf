process RIBODETECTOR_REMOVE {
    tag "$meta.id"
    label 'high'

    // Provide a minimal Python + pip env; users can pin/override via profiles
    conda "bioconda::ribodetector"

    publishDir path: "${params.outdir}/ribodetector_remove", mode: 'link', saveAs: {
        filename -> if (filename.endsWith('.fq')) return "fastq/$filename"
        else if (filename.endsWith('.log')) return "logs/$filename"
        else null
    }

    input:
    tuple val(meta), path(reads)
    path rrna_ref

    output:
    tuple val(meta), path("*.fq"), emit: fastq
    tuple val(meta), path("*_ribodetector_alignment_stats.log"), emit: log
    path "versions.yml", emit: versions

    when:
    true

    script:
    def args = task.ext.args ?: (params.ribodetector_args ?: '')
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ribodetector -t 20 \
    -l 50 \
    -i $reads \
    -m 10 \
    -e rrna \
    --chunk_size 256 \
    -o ${prefix}_norrna.fq \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector: \$(command -v ribodetector >/dev/null 2>&1 && ribodetector --version 2>/dev/null | head -n1 || echo "not_available")
        python: \$(python --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_lessrRNA.fq
    echo "stub run" > ${prefix}_ribodetector_alignment_stats.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector: stub
        python: stub
    END_VERSIONS
    """
}

