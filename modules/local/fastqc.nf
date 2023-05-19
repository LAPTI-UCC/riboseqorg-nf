
// Example module that runs fastqc on a single fastq file.
// Results are publushed to the study directory in the fastqc directory

process FASTQC {

    tag 'medium'

	publishDir "${params.output_dir}/fastqc", mode: 'copy'
	
	input:
	    file fastq 

	output:
	    path "*_fastqc.{zip,html}", emit: fastqc_full_reports

    script:
        """
        fastqc -q $fastq 
        """
}