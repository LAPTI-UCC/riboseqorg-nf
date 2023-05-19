
process FASTQ_DL {
    tag 'medium'

    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*.fastq.gz'

    input:
        val(accession)

    output:
        file "*.fastq.gz"

    script:
        """
        fastq-dl -a $accession --cpus 5
        """
}