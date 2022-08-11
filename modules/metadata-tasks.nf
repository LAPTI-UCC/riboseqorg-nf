
process GET_GSE_REPORT {

    errorStrategy {sleep 1 : 'retry'}
        
    input:
    val GSE_WNL

    output:
    path "*.xml"

	script: 
/// slicing the GSE so it does not have the /n inside (the /n is added by the splitText operator, see workflow)
    GSE = "${GSE_WNL[0..-2]}"
/// sleep_GSE introduces a random delay in the download of the files.
    sleep_GSE = "${GSE[-1]}"
    if (sleep_GSE == "0" || sleep_GSE == "1" || sleep_GSE == "3" ){
        sleep_GSE = "5"
    }

	"""
    sleep ${sleep_GSE}
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/${GSE[0..-4]}nnn/${GSE}/miniml/${GSE}_family.xml.tgz
    tar -zxvf ${GSE}_family.xml.tgz
    """
}

process GET_CSV_FROM_XML {
    /// VERY provisional, I just need to  check the result (and thus having a single folder makes it easier)
    /// The publishing directory will be the one of the study in exam
    publishDir "/home/121109636/CSV_reports"

    input:
    path xml_report

    output:
    path "*.csv"

    script:
    """
    python3 $projectDir/../scripts/xml_parsing.py  ${xml_report}
    """
}


workflow {

    params.path_to_txt = 
    input = Channel
        .fromPath("/home/121109636/CSV_reports/GSEs.txt")
        .splitText()

    GET_GSE_REPORT          ( input )
    GET_CSV_FROM_XML        ( GET_GSE_REPORT.out )
}
