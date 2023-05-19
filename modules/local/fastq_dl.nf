
process FASTQ_DL {
    tag 'medium'

    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*.fastq.gz'

    input:
        tuple val(accession), val(GSE)

    output:
        file "*.fastq.gz"

    script:
        """
        fastq-dl -a $accession --cpus 5
        """
}