FROM condaforge/mambaforge:latest

LABEL description="BAM to SQLite converter with pysam and sqlitedict" \
      author="Jack Curragh" \
      version="1.0"

# Install required packages
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.11 \
    pysam=0.23.3 \
    sqlite \
    sqlitedict \
 && mamba clean -a -y

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
