FROM condaforge/mambaforge:latest

LABEL description="Python 3.9 with pandas, sqlitedict, and numpy for data processing" \
      author="Jack Curragh" \
      version="1.0"

# Install required packages
RUN mamba install -y -c conda-forge \
    python=3.9 \
    pandas \
    numpy \
    sqlite \
    sqlitedict \
 && mamba clean -a -y

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
