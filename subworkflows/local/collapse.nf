


// Include the necessary modules for this workflow
include { FASTQC } from '../../modules/local/fastqc/main.nf'
include { FASTQ_DL } from '../../modules/local/fastq_dl/main.nf'
include { MULTIQC } from '../../modules/local/multiqc/main.nf'
include { FASTP } from '../../modules/local/fastp/main.nf'
include { COLLAPSE_FASTQ } from '../../modules/local/collapse/main.nf'
include { COLLAPSE_FASTQ as COLLAPSE_FASTQ_RAW } from '../../modules/local/collapse/main.nf'
include { COLLAPSE_FASTQ as COLLAPSE_FASTQ_TRIMMED } from '../../modules/local/collapse/main.nf'
include { FIND_ADAPTERS } from '../../modules/local/find_adapters/main.nf'
include { DETECT_ARCHITECTURE } from '../../modules/local/getRPF/detect_architecture/main.nf'

workflow collapse {
    take:
        sample_ch // Channel: [ meta, run_accession ]

    main:

        // Define the path to the adapter list
        adapter_list = file("$projectDir/resources/adapter_list.tsv")

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

        // Run FastQC on downloaded FastQ files (always run for QC purposes)
        FASTQC(fastq_files, params.adapter_list)

        // Branch workflow based on architecture detection setting
        use_arch_detect = params.use_architecture_detection ?: false

        if (use_arch_detect) {
            // Architecture detection workflow
            // Step 1: Collapse raw reads first (needed for architecture detection)
            COLLAPSE_FASTQ_RAW(FASTQ_DL.out.fastq)

            // Step 2: Detect architecture on collapsed reads
            DETECT_ARCHITECTURE(COLLAPSE_FASTQ_RAW.out.collapsed_fastq)
            
            fastq_with_adapters = FASTQ_DL.out.fastq
                .join(DETECT_ARCHITECTURE.out.adapters, by: 0)
            // Step 4: Use FASTP with seqspec-derived adapters (same as traditional workflow!)
            FASTP(fastq_with_adapters)

            // Step 5: Collapse the trimmed reads for final output
            COLLAPSE_FASTQ_TRIMMED(FASTP.out.trimmed_fastq)
            final_collapsed = COLLAPSE_FASTQ_TRIMMED.out.collapsed_fastq

            multiqc_files = Channel.empty()
                .mix(FASTQC.out.zip.map { meta, file -> file })
                .mix(FASTP.out.json_provided.map { meta, file -> file })

            versions_ch = Channel.empty()
                .mix(FASTQ_DL.out.versions)
                .mix(FASTQC.out.versions)
                .mix(DETECT_ARCHITECTURE.out.versions)
                .mix(FASTP.out.versions)
        } else {
            // Traditional adapter finding workflow
            FIND_ADAPTERS(FASTQ_DL.out.fastq, FASTQC.out.txt)

            fastq_with_adapters = FASTQ_DL.out.fastq
                .join(FIND_ADAPTERS.out.adapter_report, by: 0)

            // Run FastP for trimming
            FASTP(fastq_with_adapters)

            // Collapse FastQ files
            COLLAPSE_FASTQ(FASTP.out.trimmed_fastq)
            final_collapsed = COLLAPSE_FASTQ.out.collapsed_fastq
            
            multiqc_files = Channel.empty()
                .mix(FASTQC.out.zip.map { meta, file -> file })
                .mix(FASTP.out.json_provided.map { meta, file -> file })

            versions_ch = Channel.empty()
                .mix(FASTQ_DL.out.versions)
                .mix(FASTQC.out.versions)
                .mix(FASTP.out.versions)
        }

        // Run MultiQC only if there are input files
        // MULTIQC(multiqc_files.collect())

    emit:
        collapsed_fastq = final_collapsed
        fastqc_results  = FASTQC.out.html
        // multiqc_report  = MULTIQC.out.report
        versions        = versions_ch.collect()
}