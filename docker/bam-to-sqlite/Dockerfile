FROM condaforge/mambaforge:23.11.0-0  # pin this

LABEL description="BAM to SQLite converter with pysam and sqlitedict" \
      author="Jack Curragh" \
      version="1.0"

# Create isolated env
RUN mamba create -y -n app \
    -c conda-forge -c bioconda \
    python=3.11 \
    pysam=0.23.3 \
    sqlite \
    sqlitedict \
 && mamba clean -a -y

ENV PATH=/opt/conda/envs/app/bin:$PATH

WORKDIR /data
CMD ["/bin/bash"]
