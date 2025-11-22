FROM condaforge/mambaforge:latest

LABEL description="BAM to BED converter with pysam, samtools, numpy, and biopython" \
      author="Jack Curragh" \
      version="1.0"

# Install required packages
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.11 \
    pysam=0.23.3 \
    samtools=1.20 \
    numpy \
    biopython \
 && mamba clean -a -y

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
