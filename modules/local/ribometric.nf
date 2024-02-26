process RIBOMETRIC {

    tag 'medium'

	publishDir "${params.study_dir}/ribometric", mode: 'copy'

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

	
	input:
	    file transcriptome_bam
        file transcriptome_bam_bai

	output:
	    path "*RiboMetric.html", emit: ribometric_reports
        path "*RiboMetric.json", emit: ribometric_json

    script:
        """
        ribometric run --bam ${transcriptome_bam} --annotation ${params.annotation} --threads $task.cpus --html --json --csv > /dev/null
        """
}
