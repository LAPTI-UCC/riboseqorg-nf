FROM condaforge/mambaforge:latest

LABEL description="RiboSeq-DP-Tools container for data processing" \
      author="Jack Curragh" \
      version="1.0"

# Install conda packages (prebuilt binaries, avoids compilation)
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10.13 \
    pip \
    biopython \
    psutil \
    pysam \
 && mamba clean -a -y

# Install RiboSeq-DP-Tools from PyPI (use --no-deps, avoid pip overriding conda libs)
RUN pip install --no-cache-dir --no-deps RiboSeq-DP-Tools==0.1.10

WORKDIR /data

CMD ["/bin/bash"]
