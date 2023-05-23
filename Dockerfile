FROM ubuntu:latest

# Install necessary packages
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    default-jre \
    perl

# Download and install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    ln -s /FastQC/fastqc /usr/local/bin/fastqc

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Set PATH for Conda
ENV PATH="/opt/miniconda/bin:$PATH"

# Create a new Conda environment
RUN conda update -n base -c defaults conda
RUN conda create -n myenv

# Activate the Conda environment
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Install fastq-dl from Bioconda (specify explicit version)
RUN conda install -c bioconda fastq-dl=2.0.1

# Install SRA Toolkit and fasterq-dump
RUN apt-get install -y sra-toolkit

# Set the default command
CMD ["ls"]
