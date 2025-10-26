

process FILTER_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    publishDir "${params.outdir}/filtered_bams", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.unique_juntion_mapped.bam"),       emit: unique_juntion_mapped
    tuple val(meta), path("${prefix}.unique_non_juntion_mapped.bam"),   emit: unique_non_juntion_mapped
    tuple val(meta), path("${prefix}.max_*mm_juntion_mapped.bam"),      emit: max_mm_juntion_mapped
    tuple val(meta), path("${prefix}.max_*mm_non_juntion_mapped.bam"),  emit: max_mm_non_juntion_mapped
    tuple val(meta), path("${prefix}*.bam"),                            emit: all_bams
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    max_mm = params.max_multimappers

    // Build awk pattern for max multimappers dynamically
    // For numbers 1-9, use [1-9], for 10+, enumerate each two-digit number
    def max_mm_pattern = ""
    if (max_mm <= 9) {
        max_mm_pattern = "(\\\$0 ~ /^@/) || (\\\$0 ~ /\\\\tNH:i:[1-${max_mm}]\\\\t/ || \\\$0 ~ /\\\\tNH:i:[1-${max_mm}]\\\$/)"
    } else {
        // For 10+, need to match both single digits [1-9] and two digits separately
        def two_digit_patterns = []
        for (int i = 10; i <= max_mm; i++) {
            two_digit_patterns << "\\\\tNH:i:${i}\\\\t"
            two_digit_patterns << "\\\\tNH:i:${i}\\\$"
        }
        def two_digit_regex = two_digit_patterns.join("/ || \\\$0 ~ /")
        max_mm_pattern = "(\\\$0 ~ /^@/) || (\\\$0 ~ /\\\\tNH:i:[1-9]\\\\t/ || \\\$0 ~ /\\\\tNH:i:[1-9]\\\$/ || \\\$0 ~ /${two_digit_regex}/)"
    }

    """
    # First create the unique and max_${max_mm}mm mappers BAMs
    # Unique mappers (exactly 1)
    samtools view -h ${bam} | awk '(\$0 ~ /^@/) || (\$0 ~ /\\tNH:i:1\\t/ || \$0 ~ /\\tNH:i:1\$/)' | samtools view -bS - > temp_unique.bam
    # Max ${max_mm} mappers (NH:i:1 through NH:i:${max_mm})
    samtools view -h ${bam} | awk '${max_mm_pattern}' | samtools view -bS - > temp_max_mm.bam
    
    # Then filter each for juntion_mapped/non_juntion_mapped
    # Junction mapped: keep reads with N in CIGAR (spliced)
    samtools view -h temp_unique.bam | awk '(\$0 ~ /^@/) || (\$0 !~ /^@/ && \$6 ~ /N/)' | samtools view -bS - > ${prefix}.unique_juntion_mapped.bam
    samtools view -h temp_unique.bam | awk '(\$0 ~ /^@/) || (\$0 !~ /^@/ && \$6 !~ /N/)' | samtools view -bS - > ${prefix}.unique_non_juntion_mapped.bam

    samtools view -h temp_max_mm.bam | awk '(\$0 ~ /^@/) || (\$0 !~ /^@/ && \$6 ~ /N/)' | samtools view -bS - > ${prefix}.max_${max_mm}mm_juntion_mapped.bam
    samtools view -h temp_max_mm.bam | awk '(\$0 ~ /^@/) || (\$0 !~ /^@/ && \$6 !~ /N/)' | samtools view -bS - > ${prefix}.max_${max_mm}mm_non_juntion_mapped.bam

    # Clean up temp files
    rm temp_unique.bam temp_max_mm.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    max_mm = params.max_multimappers
    """
    touch ${prefix}.unique_juntion_mapped.bam
    touch ${prefix}.unique_non_juntion_mapped.bam
    touch ${prefix}.max_${max_mm}mm_juntion_mapped.bam
    touch ${prefix}.max_${max_mm}mm_non_juntion_mapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}