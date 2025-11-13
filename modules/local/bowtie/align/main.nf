


process BOWTIE_ALIGN_SORT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bowtie=1.3.1 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:c84c7c55c45af231883d9ff4fe706ac44c479c36-0' :
        'biocontainers/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:c84c7c55c45af231883d9ff4fe706ac44c479c36-0' }"

    // Add publishing directives
    publishDir path: "${params.outdir}/bowtie_align", mode: 'link', saveAs: { 
        filename -> if (filename.endsWith('.bam')) return "bam/$filename"
        else if (filename.endsWith('.log')) return "logs/$filename" 
        else null 
    }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam
    tuple val(meta), path("*_bowtie_alignment_stats.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readInput = meta.single_end ? "-q ${reads[0]}" : "-q -1 ${reads[0]} -2 ${reads[1]}"

    """
    INDEX=\$(find -L ./ -name "*.3.ebwt" | sed 's/\\.3.ebwt\$//')

    bowtie \
        -S -a --norc \
        -p ${task.cpus} \
        -v 3 \
        --seedlen 25 \
        -x \$INDEX \
        -f ${reads} 2> ${meta.id}_bowtie_alignment_stats.log | samtools view -bS - | samtools sort -@ ${task.cpus} -o ${prefix}.bam


    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${meta.id}_bowtie_alignment_stats.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: 1.3.1
        samtools: 1.21
    END_VERSIONS
    """
}