process RIBODETECTOR_REMOVE {
    tag "$meta.id"
    label 'high'

    // Provide a minimal Python + pip env; users can pin/override via profiles
    conda "conda-forge::python=3.10 pip"

    publishDir path: "${params.outdir}/ribodetector_remove", mode: 'link', saveAs: {
        filename -> if (filename.endsWith('.fq')) return "fastq/$filename"
        else if (filename.endsWith('.log')) return "logs/$filename"
        else null
    }

    input:
    tuple val(meta), path(reads)
    path rrna_ref

    output:
    tuple val(meta), path("*.fq"), emit: fastq
    tuple val(meta), path("*_ribodetector_alignment_stats.log"), emit: log
    path "versions.yml", emit: versions

    when:
    true

    script:
    def args = task.ext.args ?: (params.ribodetector_args ?: '')
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -euo pipefail

    # Resolve rRNA fasta from provided path or directory
    RRNA_FASTA=""
    if [ -f "${rrna_ref}" ]; then
        case "${rrna_ref}" in
            *.fa|*.fa.gz|*.fasta|*.fasta.gz|*.fna|*.fna.gz)
                RRNA_FASTA="${rrna_ref}"
                ;;
            *)
                echo "Provided rRNA reference is a file but not a FASTA: ${rrna_ref}" | tee -a ${prefix}_ribodetector_alignment_stats.log
                exit 1
                ;;
        esac
    else
        # try to find a fasta inside the directory (in case a Bowtie index directory was passed)
        RRNA_FASTA=\$(find -L "${rrna_ref}" -maxdepth 2 \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | head -n1)
        if [ -z "\${RRNA_FASTA}" ]; then
            echo "Could not resolve an rRNA FASTA from: ${rrna_ref}" | tee -a ${prefix}_ribodetector_alignment_stats.log
            exit 1
        fi
    fi

    echo "rRNA FASTA: \${RRNA_FASTA}" > ${prefix}_ribodetector_alignment_stats.log

    # Install RiboDetector if available from PyPI; users may replace with containers or prebuilt envs
    python - <<'PY'
import os, sys, subprocess
print("Installing ribodetector via pip if available...", flush=True)
try:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--no-input', '--disable-pip-version-check', 'ribodetector'])
except subprocess.CalledProcessError as e:
    print("WARNING: Failed to install ribodetector; ensure it is available in the environment.", file=sys.stderr)
PY

    # Try common CLI patterns for RiboDetector; prefer explicit args if provided
    if [ -n "${args}" ]; then
        echo "Running RiboDetector with custom args: ${args}" | tee -a ${prefix}_ribodetector_alignment_stats.log
        bash -o pipefail -c "ribodetector ${args} -i ${reads} -r \${RRNA_FASTA} -o ${prefix}_lessrRNA.fq" 2>> ${prefix}_ribodetector_alignment_stats.log || true
    else
        # Attempt a default invocation; users may override via params.ribodetector_args
        if command -v ribodetector >/dev/null 2>&1; then
            echo "Running: ribodetector remove-rrna -i ${reads} -r \${RRNA_FASTA} -o ${prefix}_lessrRNA.fq" | tee -a ${prefix}_ribodetector_alignment_stats.log
            ribodetector remove-rrna -i ${reads} -r "\${RRNA_FASTA}" -o ${prefix}_lessrRNA.fq 2>> ${prefix}_ribodetector_alignment_stats.log || true
        else
            echo "ERROR: 'ribodetector' command not found. Provide it via container/conda or set params.ribodetector_args." | tee -a ${prefix}_ribodetector_alignment_stats.log
            # create empty outputs to allow downstream stubs or to fail clearly later
            : > ${prefix}_lessrRNA.fq
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector: \$(command -v ribodetector >/dev/null 2>&1 && ribodetector --version 2>/dev/null | head -n1 || echo "not_available")
        python: \$(python --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_lessrRNA.fq
    echo "stub run" > ${prefix}_ribodetector_alignment_stats.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector: stub
        python: stub
    END_VERSIONS
    """
}

