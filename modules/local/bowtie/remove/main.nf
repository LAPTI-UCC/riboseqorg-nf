

process BOWTIE_REMOVE {
    tag "$meta.id"
    label 'high'

    conda "bioconda::bowtie=1.3.1 bioconda::samtools=1.21"

    // Add publishing directives
    publishDir path: "${params.outdir}/bowtie_remove", mode: 'link', saveAs: { 
        filename -> if (filename.endsWith('.fq')) return "fastq/$filename"
        else if (filename.endsWith('.log')) return "logs/$filename" 
        else null 
    }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.fq"), emit: fastq
    tuple val(meta), path("*_bowtie_alignment_stats.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readInput = meta.single_end ? "-q ${reads[0]}" : "-q -1 ${reads[0]} -2 ${reads[1]}"

    """
    INDEX=\$(find -L ./ -name "*.3.ebwt" | sed 's/\\.3.ebwt\$//')

    bowtie \
        -v 3\
        -p 8\
        -x \$INDEX \
        -f ${reads} \
        --un ${prefix}_lessrRNA.fq 1> /dev/null  2>  ${prefix}_bowtie_alignment_stats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_lessrRNA.fq
    touch ${prefix}_bowtie_alignment_stats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: 1.3.1
    END_VERSIONS
    """
}