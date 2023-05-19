process FIND_ADAPTERS {
    publishDir "$projectDir/$params.data_dir/$params.GSE/fastq", mode: 'copy', pattern: '*_adpater_report.tsv'

    input:
        file raw_fastq

    output:
        file "${raw_fastq}_adpater_report.tsv"

    script:
        """
        python3 $projectDir/scripts/get_adapters.py -q $raw_fastq -o "${raw_fastq}_adpater_report.tsv"
        """
}