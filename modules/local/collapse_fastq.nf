

process COLLAPSE_FASTQ {
	publishDir "${params.study_dir}/collapsed_fastq", mode: 'copy'

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        path fastq

    output:
        path "*.fastq.gz"

    script:
    """
    python3 $projectDir/scripts/collapse_fastq.py -i $fastq -o ${fastq.baseName}_collapsed.fastq.gz
    """
}