process RIBOMETRIC {

    tag 'medium'

	publishDir "${params.study_dir}/ribometric", mode: 'copy'
	
	input:
	    file transcriptome_bam
        file transcriptome_bam_bai

	output:
	    path "*RiboMetric.csv", emit: ribometric_reports
        path "*RiboMetric.json", emit: ribometric_json

    script:
        """
        RiboMetric run --bam ${transcriptome_bam} --annotation ${params.annotation} --threads $task.cpus --html --json --csv 2>&1
        """
}
