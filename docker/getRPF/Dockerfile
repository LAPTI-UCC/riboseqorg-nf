FROM condaforge/mambaforge:latest

LABEL description="getRPF container for riboseq processing" \
      author="Jack Curragh" \
      version="1.0"

# Install conda packages (includes pre-compiled binaries to avoid build issues)
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    biopython \
    git && \
    mamba clean -a -y

# Install getRPF from GitHub
RUN pip install git+https://github.com/JackCurragh/get-RPF.git

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
