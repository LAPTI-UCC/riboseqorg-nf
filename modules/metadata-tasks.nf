
process GET_GSE_REPORT {
        
    input:
    val GSE_WNL

    output:
    path "*.xml.tgz",emit: xml_tgz
    val GSE, emit: GSE

	script: 
/// slicing the GSE so it does not have the /n inside (the /n is added by the splitText operator, see workflow)
    GSE = "${GSE_WNL[0..-2]}"
/// sleep_GSE introduces a random delay in the download of the files.
    sleep_GSE = "${GSE[-1]}"
    if (sleep_GSE == "0"){
        sleep_GSE = "2"
    }

	"""
    sleep ${sleep_GSE}
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/${GSE[0..-4]}nnn/${GSE}/miniml/${GSE}_family.xml.tgz
    """
}

process EXTRACT_XML_REPORT {

    errorStrategy { task.exitStatus == 2 ? 'retry' : 'terminate' }

    input:
    file compressed_xml
    val GSE

    output:
    path "*.xml"

    script:
    """
    tar -zxvf ${compressed_xml} ${GSE}_family.xml 
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
    /*
    input = Channel.of( "GSE180669","GSE156796", "GSE185286", "GSE185458","GSE158141",
"GSE173856","GSE136940","GSE130465","GSE157361","GSE157423","GSE152554","GSE152556",
"GSE152558","GSE167704","GSE167223","GSE166874","GSE160917","GSE144539","GSE158881",
"GSE157063","GSE134152","GSE143301","GSE126660","GSE123981","GSE131074","GSE140565",
"GSE129757","GSE137626","GSE138278","GSE128538","GSE139399","GSE126736","GSE139880",
"GSE104503","GSE105782","GSE128344","GSE119681","GSE132725","GSE125725", "GSE133125",
"GSE133111","GSE112705","GSE116233","GSE127713","GSE121952","GSE110618","GSE123675",
"GSE97286","GSE102216") */

    params.path_to_txt = 
    input = Channel
        .fromPath("/home/121109636/CSV_reports/GSEs.txt")
        .splitText()

    GET_GSE_REPORT          ( input )
    EXTRACT_XML_REPORT      ( GET_GSE_REPORT.out.xml_tgz, GET_GSE_REPORT.out.GSE )
    GET_CSV_FROM_XML        ( EXTRACT_XML_REPORT.out )
}
