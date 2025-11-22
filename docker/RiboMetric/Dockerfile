FROM condaforge/mambaforge:latest

LABEL description="RiboMetric container for ribosome profiling QC" \
      author="Jack Curragh" \
      version="1.0"

# Install conda packages (precompiled, avoids compilation issues)
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    biopython \
    pysam \
    git \
 && mamba clean -a -y

# Install RiboMetric from GitHub safely
RUN pip install --no-cache-dir --no-deps \
    git+https://github.com/JackCurragh/RiboMetric.git

WORKDIR /data

CMD ["/bin/bash"]
