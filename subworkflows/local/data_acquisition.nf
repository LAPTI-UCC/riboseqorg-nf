

include { fetch_data } from './fetch_data'

workflow DATA_ACQUISITION {
    take:
    sample_sheet

    main:
    samples_ch = Channel
        .fromPath(sample_sheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            def meta = [
                id: row.Run,
                study_accession: row.study_accession,
            ]
            [ meta, row.Run ]
        }
    fetch_data(samples_ch)

    emit:
    samples = fetch_data.out
}