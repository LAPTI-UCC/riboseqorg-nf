FROM condaforge/mambaforge:latest

LABEL description="RiboMetric container for ribosome profiling QC" \
      author="Jack Curragh" \
      version="1.0"

# Install conda packages
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    biopython \
    pysam && \
    mamba clean -a -y

# Install RiboMetric from GitHub
RUN pip install git+https://github.com/JackCurragh/RiboMetric.git

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
