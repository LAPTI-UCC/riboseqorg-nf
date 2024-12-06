

include { RIBOMETRIC } from '../../modules/local/ribometric/main'

workflow ANALYSIS {
    take:
    processed_data
    ribometric_annotation

    main:
    RIBOMETRIC(processed_data, ribometric_annotation)

    emit:
    ribometric_html = RIBOMETRIC.out.html
    ribometric_csv = RIBOMETRIC.out.csv
    ribometric_json = RIBOMETRIC.out.json

}