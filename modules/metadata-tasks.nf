
process GET_GSE_REPORT {
        
    input:
    val GSE_WNL

    output:
    path "*.xml"

	script: 
/// slicing the GSE so it does not have the /n inside (the /n is added by the splitText operator, see workflow)
    GSE = "${GSE_WNL[0..-2]}"
/// sleep ${GSE[-1]} introduces a random delay in the download of the files.
	"""
    sleep ${GSE[-1]}
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/${GSE[0..-4]}nnn/${GSE}/miniml/${GSE}_family.xml.tgz
    tar -xzvf ${GSE}_family.xml.tgz
    """
}

process GET_CSV_FROM_XML {
    /// VERY provisional, I just need to iterate this process over several GSE and check the result
    /// The publishing directory will be the one of the study
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
        .fromPath("/home/121109636/CSV_reports/GSEs2.txt")
        .splitText()
        .view()

    GET_GSE_REPORT          ( input )
    GET_CSV_FROM_XML        ( GET_GSE_REPORT.out )
}
