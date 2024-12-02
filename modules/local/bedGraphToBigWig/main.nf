process BEDGRAPH_TO_BIGWIG {

	publishDir "$params.outdir/bigwigs", mode: 'copy', pattern: '*.bw'

	input:
	file bedfile /// from sorted_beds ///
	
    output:
	file "*.bw"  /// into bigwigs ///

	"""
	$projectDir/scripts/bedGraphToBigWig ${bedfile} $params.chrom_sizes_file ${bedfile.baseName}.coverage.bw
	"""
}
