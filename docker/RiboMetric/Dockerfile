FROM condaforge/mambaforge:latest

LABEL description="RiboMetric container for ribosome profiling QC" \
      author="Jack Curragh" \
      version="1.0"

# Install conda packages (precompiled, avoids compilation issues)
# Install pyarrow from conda-forge to get proper shared libraries
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    biopython \
    pysam \
    git \
    rich \
    pyarrow \
    samtools \
 && mamba clean -a -y

# Install RiboMetric from GitHub safely
# Use --no-deps for packages already installed via conda to avoid conflicts
RUN pip install --no-build-isolation \
    git+https://github.com/JackCurragh/RiboMetric.git

WORKDIR /data

CMD ["/bin/bash"]
