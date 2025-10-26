include { CHECK_CLEANLINESS } from '../../modules/local/getRPF/check_cleanliness/main'

workflow QUALITY_CONTROL {
    take:
    samples  // channel: [meta, input_file]

    main:
    ch_versions = Channel.empty()

    // Cleanliness check
    CHECK_CLEANLINESS(samples, params.collapsed_read_header_pattern)
    ch_versions = ch_versions.mix(CHECK_CLEANLINESS.out.versions)

    emit:
    clean_samples = samples                                             // channel: [meta, file]
    versions      = ch_versions                                         // channel: [versions.yml]
}