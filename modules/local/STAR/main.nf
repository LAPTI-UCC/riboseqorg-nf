

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
    tuple val(meta), path("*.Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam, optional: true
    tuple val(meta), path("*.Log.final.out"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trim_front = params.trim_front > 0 ? "--clip3pNbases ${params.trim_front}" : ''
    def alignment_type = params.alignment_type == 'Local' ? '--alignEndsType Local' : ''
    def allow_introns = params.allow_introns ? '--alignIntronMax 1000000 --alignMatesGapMax 1000000' : ''
    def unzip_command = reads.name.endsWith('.gz') ? 'zcat' : 'cat'
    def output_transcriptome_bam = params.save_star_transcriptome_bam ? "--quantMode TranscriptomeSAM" : ""

    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${prefix}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFilterMultimapNmax ${params.max_multimap} \\
        --outFilterMismatchNmax ${params.mismatches} \\
        $output_transcriptome_bam \\
        $alignment_type \\
        $allow_introns \\
        $trim_front \\
        --readFilesCommand $unzip_command \\
        --sjdbGTFfile $gtf \\
        $args
        STAR \
     --genomeDir /path/to/star_index \          # Directory containing STAR genome index \
     --readFilesIn sample_R1.fastq.gz \         # Input fastq file \
     --readFilesCommand zcat \                  # Command to unzip input files on the fly \
     --runThreadN 8 \                           # Number of threads to use \
     --alignEndsType EndToEnd \                 # Force end-to-end read mapping, no soft-clipping \
     --outFilterMismatchNmax 2 \                # Maximum number of mismatches per read \
     --outFilterMultimapNmax 20 \               # Report up to 20 multimappings per read \
     --winAnchorMultimapNmax 50 \               # Maximum number of loci anchors are allowed to map to \
     --seedSearchStartLmax 30 \                 # Maximum length of the seed region \
     --alignSJoverhangMin 4 \                   # Minimum overhang for splice junctions \
     --alignSJDBoverhangMin 1 \                 # Minimum overhang for annotated splice junctions \
     --alignIntronMin 20 \                      # Minimum intron size \
     --alignIntronMax 100000 \                  # Maximum intron size \
     --outSAMtype BAM SortedByCoordinate \      # Output sorted BAM file \
     --outFileNamePrefix sample_ \              # Prefix for output files \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    
}

