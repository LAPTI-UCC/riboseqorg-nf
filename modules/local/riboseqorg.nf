

process BAM_TO_SQLITE {

	publishDir "$params.study_dir/sqlites", mode: 'copy', pattern: '*.sqlite'

	input:
	file sorted_bam /// from sorted_bams ///

	output:
	file "*.sqlite" /// into sqlite_ch ///

	"""
	python3 $projectDir/scripts/bam_to_sqlite.py --bam ${sorted_bam} --annotation $params.annotation_sqlite --output ${sorted_bam.baseName}.sqlite
	"""
}