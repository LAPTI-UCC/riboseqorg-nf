process RIBOMETRIC {

    tag 'medium'

	publishDir "${params.study_dir}/ribometric", mode: 'copy'
	
	input:
	    file transcriptome_bam

	output:
	    path "*_ribometric.html", emit: ribometric_reports
        path "*_ribometric.json", emit: ribometric_json

    script:
        """
        ribometric run --bam ${transcriptome_bam} --annotation ${params.annotation} --threads 5 --html --csv
        """
}
