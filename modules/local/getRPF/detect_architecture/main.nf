  process DETECT_ARCHITECTURE {
      tag "${meta.id}"
      label 'process_medium'

      conda "${projectDir}/conda/getRPF.yml"

      publishDir "${params.outdir}/getRPF/detect", mode: 'copy'

    //   errorStrategy 'ignore'

      input:
      tuple val(meta), path(input_file)

      output:
      tuple val(meta), path("*.seqspec.yaml")           , emit: seqspec
      tuple val(meta), path("*.extraction_report.json") , emit: report
      tuple val(meta), path("*_rpfs.fastq")             , emit: rpfs, optional: true
      path "versions.yml"                               , emit: versions

      script:
      def args = task.ext.args ?: ''
      def prefix = task.ext.prefix ?: "${meta.id}"
      """
      getRPF detect-architecture \\
          ${input_file} \\
          ${prefix}_rpfs.fastq \\
          -f collapsed \\
          --generate-seqspec \\
          --max-reads 5000 \\
          $args

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          getrpf: \$(getRPF --version 2>&1 | head -n1 | sed 's/^.*version //g')
      END_VERSIONS
      """

      stub:
      def prefix = task.ext.prefix ?: "${meta.id}"
      """
      touch ${prefix}_rpfs.seqspec.yaml
      touch ${prefix}_rpfs.extraction_report.json
      touch ${prefix}_rpfs.fastq

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          getrpf: 1.3.0
      END_VERSIONS
      """
  }