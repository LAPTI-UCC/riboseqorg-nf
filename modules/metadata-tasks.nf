
process GET_GSE_REPORT {

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
     
    input:
        tuple val(GSE_WNL), val(srp)

    output:
        path "*.xml"

	script: 
/// slicing the GSE so it does not have the /n inside (the /n is added by the splitText operator, see workflow)
        GSE = "${GSE_WNL[0..-1]}"
    /// sleep_GSE introduces a random delay in the download of the files.
        def z = ["4", "5", "6", "7", "8", "9"]
        Random rnd = new Random()

        Sleep_time = (z[rnd.nextInt(z.size)])

        """
        sleep ${Sleep_time}
        wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/${GSE[0..-4]}nnn/${GSE}/miniml/${GSE}_family.xml.tgz
        tar -zxvf ${GSE}_family.xml.tgz
        """
}


process GET_CSV_FROM_XML {
    // errorStrategy 'ignore'

    input:
    path xml_report

    output:
    path "*.csv"

    script:
    """
    python3 $projectDir/scripts/xml_parsing.py  ${xml_report}
    """
}


process ASSESS_LIBRARY_STRATEGY {
    publishDir "$projectDir/$params.data_dir/$params.GSE", mode: 'copy'

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    // errorStrategy 'ignore'


    input:
    path csv

    output:
    path csv

    script:
    """
    python3 $projectDir/scripts/library_strategy_definer.py ${csv}
    """

}


process CHECK_METADATA_REPORT {

    input: 
    path csv

    script:
    """
    python3 $projectDir/scripts/metadata_completeness_checker.py ${csv}
    """

}
/*
workflow {

    ///params.path_to_txt = 
    input = Channel
        .fromPath("/home/121109636/CSV_reports/GSEs.txt")
        ///.fromPath("/home/gionmattia/Desktop/ResearchProject/CSV_reports/GSEs.txt")
        .splitText()

    GET_GSE_REPORT          ( input )
    GET_CSV_FROM_XML        ( GET_GSE_REPORT.out )
    ASSESS_LIBRARY_STRATEGY ( GET_CSV_FROM_XML.out )
}
*/