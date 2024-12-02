

process STAR_ALIGN {
    tag "$meta.id"
    label 'high'

    conda "bioconda::star=2.7.10a"

    // Add publishing directives
    publishDir path: "${params.outdir}/star_align", mode: 'link', saveAs: { 
        filename -> if (filename.endsWith('toTranscriptome.out.bam')) return "transcriptome_bam/$filename" 
        else if  (filename.endsWith('.bam')) return "bam/$filename"
        else if (filename.endsWith('.out')) return "logs/$filename" else null 
        }


    input:
    tuple val(meta), path(reads)
    path index
    path gtf

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("*.Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam
    tuple val(meta), path("*.Log.final.out"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trim_front = params.trim_front > 0 ? "--clip3pNbases ${params.trim_front}" : ''
    def alignment_type = params.alignment_type == 'Local' ? '--alignEndsType Local' : ''
    def allow_introns = params.allow_introns ? '--alignIntronMax 1000000 --alignMatesGapMax 1000000' : ''
    def unzip_command = reads.name.endsWith('.gz') ? 'zcat' : 'cat'
    
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${prefix}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFilterMultimapNmax ${params.max_multimap} \\
        --outFilterMismatchNmax ${params.mismatches} \\
        --quantMode TranscriptomeSAM \\
        $alignment_type \\
        $allow_introns \\
        $trim_front \\
        --readFilesCommand $unzip_command \\
        --sjdbGTFfile $gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}