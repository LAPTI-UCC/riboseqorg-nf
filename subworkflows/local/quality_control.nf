

include { GET_RPF } from '../../modules/local/getRPF/main'

workflow QUALITY_CONTROL {
    take:
    samples

    main:
    GET_RPF(samples, params.collapsed_read_header_pattern)
    
    GET_RPF.out.rpf_checks
        .branch {
            pass: !it[1].text.contains("[FAIL]")
            fail: it[1].text.contains("[FAIL]")
        }
        .set { rpf_results }

    rpf_results.fail
        .map { meta, _rpf_checks -> 
            log.warn("Sample ${meta.id} failed RPF checks. Skipping alignment.")
            return meta.id
        }
        .set { failed_samples }

    passed_samples = rpf_results.pass
        .join(samples, by: [0])
        .map { meta, _rpf_checks, input_file -> [meta, input_file] }

    emit:
    passed_samples
    failed_samples
}