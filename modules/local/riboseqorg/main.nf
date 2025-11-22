

process BAM_TO_SQLITE {

	conda "bioconda::pysam=0.23.3 conda-forge::sqlite=3.39.3 conda-forge::sqlitedict=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/lapti-ucc/riboseqorg-nf/releases/download/containers-latest/bam-to-sqlite.sif' :
        'ghcr.io/lapti-ucc/riboseqorg-nf-bam-to-sqlite:latest' }"

	publishDir "$params.outdir/sqlites", mode: 'copy', pattern: '*.sqlite'

	input:
	file sorted_bam /// from sorted_bams ///

	output:
	file "*.sqlite" /// into sqlite_ch ///

	"""
	python3 $projectDir/scripts/bam_to_sqlite.py --bam ${sorted_bam} --annotation $params.annotation_sqlite --output ${sorted_bam.baseName}.sqlite
	"""
}


process GWIPS_INSERTS {

	conda "conda-forge::python=3.9 conda-forge::pandas conda-forge::sqlite"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/lapti-ucc/riboseqorg-nf/releases/download/containers-latest/python-pandas-sqlite.sif' :
        'ghcr.io/lapti-ucc/riboseqorg-nf-python-pandas-sqlite:latest' }"

	publishDir "$params.outdir/gwips_inserts", mode: 'copy', pattern: '*.sql'

	input:
	file run_metadata
	file study_metadata 
	file annotation_inventory_sqlite

	output:
	file "*.txt" /// into gwips_inserts ///

	"""
	python3 $projectDir/scripts/write_GWIPS_inserts.py -s ${study_metadata} -m ${run_metadata} --db ${annotation_inventory_sqlite} 
	"""
}

process TRIPS_INSERTS {

	conda "conda-forge::python=3.9 conda-forge::pandas conda-forge::sqlite"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/lapti-ucc/riboseqorg-nf/releases/download/containers-latest/python-pandas-sqlite.sif' :
        'ghcr.io/lapti-ucc/riboseqorg-nf-python-pandas-sqlite:latest' }"

	publishDir "$params.outdir/trips_inserts", mode: 'copy', pattern: '*.sql'

	input:
	file run_metadata
	file study_metadata 
	file trips_sqlite

	output:
	file "*.txt" /// into trips_inserts ///

	"""
	python3 $projectDir/scripts/write_TRIPS_inserts.py -s ${study_metadata} -r ${run_metadata} -m db -d ${annotation_inventory_sqlite} 
	"""

}