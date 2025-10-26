process BAM_TO_BED {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python=3.9 bioconda::pysam bioconda::samtools conda-forge::biopython"

    publishDir "$params.outdir/bedgraphs", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(offsets)

    output:
    tuple val(meta), path("*.sorted.bedgraph"), emit: bedgraph
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_stranded = args.contains('--stranded')
    """
    python3 $projectDir/bin/bam_to_bed.py \\
        --bam ${bam} \\
        --offsets ${offsets} \\
        --prefix ${prefix} \\
        $args

    if [ $is_stranded ]; then
        sort -k1,1 -k2,2n ${prefix}.forward.bedgraph > ${prefix}.forward.sorted.bedgraph
        sort -k1,1 -k2,2n ${prefix}.reverse.bedgraph > ${prefix}.reverse.sorted.bedgraph
    else
        sort -k1,1 -k2,2n ${prefix}.bedgraph > ${prefix}.sorted.bedgraph
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
        sort: \$(sort --version | sed -n 1p | sed 's/^.*sort (GNU coreutils) //g')
    END_VERSIONS


    
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_stranded = task.ext.args?.contains('--stranded') ?: false
    """
    if [ $is_stranded ]; then
        touch ${prefix}.forward.sorted.bedgraph
        touch ${prefix}.reverse.sorted.bedgraph
    else
        touch ${prefix}.sorted.bedgraph
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.5
        pysam: 0.19.0
        sort: 8.32
    END_VERSIONS
    """
}