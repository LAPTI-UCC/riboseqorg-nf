
process GET_GSE_REPORT {

    errorStrategy 'retry'
        
    input:
    val GSE_WNL

    output:
    path "*.xml"

	script: 
/// slicing the GSE so it does not have the /n inside (the /n is added by the splitText operator, see workflow)
    GSE = "${GSE_WNL[0..-2]}"
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
