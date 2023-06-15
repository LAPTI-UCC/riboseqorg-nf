

process FASTQC {

    tag 'medium'

	publishDir "${params.study_dir}/fastqc", mode: 'copy'
	
	input:
	    file fastq 

	output:
	    path "*_fastqc.html", emit: fastqc_html
        path "*/fastqc_data.txt", emit: fastqc_data

    script:
        """
        fastqc -q $fastq 
        unzip *_fastqc.zip
        """
}