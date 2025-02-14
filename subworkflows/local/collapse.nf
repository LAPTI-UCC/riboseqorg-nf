


// Include the necessary modules for this workflow
include { FASTQC } from '../../modules/local/fastqc/main.nf'
include { FASTQ_DL } from '../../modules/local/fastq_dl/main.nf'
include { MULTIQC } from '../../modules/local/multiqc/main.nf'
include { FASTP } from '../../modules/local/fastp/main.nf'
include { COLLAPSE_FASTQ } from '../../modules/local/collapse/main.nf'
include { FIND_ADAPTERS } from '../../modules/local/find_adapters/main.nf'

workflow collapse {
    take:
        sample_ch // Channel: [ meta, run_accession ]

    main:

        // Define the path to the adapter list
        adapter_list = file("$projectDir/scripts/adapter_list.tsv")

        // Download FastQ files
        FASTQ_DL(sample_ch)
        fastq_files = FASTQ_DL.out.fastq
            .flatMap { meta, files ->
                if (files instanceof Collection) {
                    files.collect { file -> tuple(meta, file) }
                } else {
                    [tuple(meta, files)]
                }
            }

        // Run FastQC on downloaded FastQ files
        FASTQC(fastq_files, adapter_list)

        FIND_ADAPTERS(FASTQ_DL.out.fastq, FASTQC.out.txt)
        // Run FastP for trimming
        FASTP(FASTQ_DL.out.fastq, FIND_ADAPTERS.out.adapter_report)

        multiqc_files = Channel.empty()
        multiqc_files = multiqc_files.mix(FASTQC.out.zip.map { meta, file -> file })
        multiqc_files = multiqc_files.mix(FASTP.out.json_provided.map { meta, file -> file })
    // Add other QC outputs as needed

        // Run MultiQC only if there are input files
        // MULTIQC(multiqc_files.collect())

        // Collapse FastQ files
        COLLAPSE_FASTQ(FASTP.out.trimmed_fastq)

    emit:
        collapsed_fastq = COLLAPSE_FASTQ.out.collapsed_fastq
        fastqc_results  = FASTQC.out.html
        fastp_results   = FASTP.out.json_provided
        // multiqc_report  = MULTIQC.out.report
        versions        = Channel.empty()
                            .mix(FASTQ_DL.out.versions)
                            .mix(FASTQC.out.versions)
                            .mix(FASTP.out.versions)
                            // .mix(MULTIQC.out.versions)
                            .collect()
}