
process GET_GSE_REPORT {
        
    input:
    val GSE

    output:
    path "*.xml"

	script: 
/// sleep ${GSE[-1]} introduces a random delay in the download of the files.
	"""
    sleep ${GSE[-1]} 
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/${GSE[0..-4]}nnn/${GSE}/miniml/${GSE}_family.xml.tgz
    tar -xzvf ${GSE}_family.xml.tgz
    """
}

/*
workflow {

    input = Channel.of( "GSE180669" ) 
    GET_GSE_REPORT          ( input )
    GET_GSE_REPORT.out.view()   
}
*/