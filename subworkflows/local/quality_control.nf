include { CHECK_CLEANLINESS } from '../../modules/local/getRPF/check_cleanliness/main'
include { DETECT_ARCHITECTURE } from '../../modules/local/getRPF/detect_architecture/main'

workflow QUALITY_CONTROL {
    take:
    samples  // channel: [meta, input_file]

    main:
    ch_versions = Channel.empty()

    // Initial cleanliness check
    CHECK_CLEANLINESS(samples, params.collapsed_read_header_pattern)
    ch_versions = ch_versions.mix(CHECK_CLEANLINESS.out.versions)

    // Branch samples based on initial check results
    CHECK_CLEANLINESS.out.rpf_checks
        .join(samples, by: [0])
        .branch { meta, rpf_checks, input_file ->
            pass: !rpf_checks.text.contains("[FAIL]")
                return [meta, input_file]
            fail: rpf_checks.text.contains("[FAIL]")
                // log.warn("Sample ${meta.id} failed initial RPF checks. Attempting architecture detection.")
                return [meta, input_file]
        }
        .set { initial_results }

    // Process failed samples through architecture detection to get seqspec
    DETECT_ARCHITECTURE(initial_results.fail)
    ch_versions = ch_versions.mix(DETECT_ARCHITECTURE.out.versions)

    emit:
    clean_samples           = initial_results.pass                      // channel: [meta, file] - samples that passed initial check
    failed_samples          = initial_results.fail                      // channel: [meta, file] - samples that failed initial check
    discovered_seqspecs     = DETECT_ARCHITECTURE.out.seqspec           // channel: [meta, seqspec] - seqspecs for failed samples
    architecture_reports    = DETECT_ARCHITECTURE.out.report            // channel: [meta, report] - detection reports
    versions                = ch_versions                               // channel: [versions.yml]
}