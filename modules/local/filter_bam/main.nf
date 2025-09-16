

process FILTER_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    publishDir "${params.outdir}/filtered_bams", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.unique_juntion_mapped.bam"),       emit: unique_juntion_mapped
    tuple val(meta), path("${prefix}.unique_non_juntion_mapped.bam"),   emit: unique_non_juntion_mapped
    tuple val(meta), path("${prefix}.max_10mm_juntion_mapped.bam"),          emit: max_10mm_juntion_mapped
    tuple val(meta), path("${prefix}.max_10mm_non_juntion_mapped.bam"),      emit: max_10mm_non_juntion_mapped
    tuple val(meta), path("${prefix}*.bam"),                 emit: all_bams
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # First create the unique and max_10mm mappers BAMs
    samtools view -h ${bam} | awk '(\$0 ~ /^@/) || (\$0 ~ /NH:i:1/)' | samtools view -bS - > temp_unique.bam
    samtools view -h ${bam} | awk '(\$0 ~ /^@/) || (\$0 ~ /NH:i:1/ || (\$0 ~ /HI:i:1/))' | samtools view -bS - > temp_max_10mm.bam
    
    # Then filter each for juntion_mapped/non_juntion_mapped
    samtools view -h temp_unique.bam -o ${prefix}.unique_juntion_mapped.bam
    samtools view -h temp_unique.bam | grep -v -P "\\t[0-9]+M[0-9]+N" | samtools view -bS - > ${prefix}.unique_non_juntion_mapped.bam
    
    samtools view -h temp_max_10mm.bam -o ${prefix}.max_10mm_juntion_mapped.bam
    samtools view -h temp_max_10mm.bam | grep -v -P "\\t[0-9]+M[0-9]+N" | samtools view -bS - > ${prefix}.max_10mm_non_juntion_mapped.bam
    
    # Clean up temp files
    rm temp_unique.bam temp_max_10mm.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unique_junction_mapped.bam
    touch ${prefix}.unique_not_junction_mapped.bam
    touch ${prefix}.max_10mm_junction_mapped.bam
    touch ${prefix}.max_10mm_not_junction_mapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}