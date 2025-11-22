FROM condaforge/mambaforge:latest

LABEL description="getRPF container for riboseq processing" \
      author="Jack Curragh" \
      version="1.0"

RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    biopython \
 && mamba clean -a -y

# Install getRPF (pure Python → safe with pip)
RUN pip install --no-cache-dir git+https://github.com/JackCurragh/get-RPF.git

WORKDIR /data

CMD ["/bin/bash"]
