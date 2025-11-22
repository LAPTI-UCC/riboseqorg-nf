FROM condaforge/mambaforge:latest

LABEL description="RiboSeq-DP-Tools container for data processing" \
      author="Jack Curragh" \
      version="1.0"

# Install conda packages including dependencies from conda to avoid compilation
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10.13 \
    pip \
    biopython \
    psutil \
    pysam && \
    mamba clean -a -y

# Install RiboSeq-DP-Tools from PyPI (dependencies already satisfied)
RUN pip install --no-deps RiboSeq-DP-Tools==0.1.10

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
