FROM condaforge/mambaforge:latest

LABEL description="RiboMetric container for ribosome profiling QC" \
      author="Jack Curragh" \
      version="1.0"

# Install conda packages (includes pre-compiled binaries to avoid build issues)
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    biopython \
    pysam \
    git && \
    mamba clean -a -y

# Install RiboMetric from GitHub
RUN pip install git+https://github.com/JackCurragh/RiboMetric.git

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
