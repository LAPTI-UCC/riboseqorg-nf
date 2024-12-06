# Ribo-Seq Data Processing
## Introduction 

This pipeline processes data to be hosted on [RiboSeq.Org](https://riboseq.org) resources. Processed files are made available through the [RiboSeq Data Portal](https://rdp.ucc.ie) and samples populate [GWIPS-viz](https://gwips.ucc.ie) and [Trips-Viz](https://tyrips.ucc.ie).

## Requirements 
This pipeline can be run using each of the following container methods
| Method      | Instructions                                                                                   |
| ----------- | ---------------------------------------------------------------------------------------------- |
| Singularity | [docs.syslabs.io](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)              |
| Docker      | [docs.docker.com](https://docs.docker.com/engine/install/)                                     |
| Conda       | [docs.conda.io](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)  |


## Setup
##### Singularity
```
sudo singularity build singularity/pipeline Singularity
```
Then as the profile `singularity` specifies `container = 'singularity/pipeline'` use the following to execute:
```
nextflow run main.nf -profile singularity
```

##### Docker
Uses mopdule specific containers.

Then as the profile `docker` specifies `container = 'pipeline-image:latest'` use the following to execute:
```
nextflow run main.nf -profile docker
```

##### Conda 

Uses module specific conda environments. 
```
nextflow run main.nf -profile conda
```

## Usage
Call the pipeline directly
```
nextflow run main.nf
```

Run with all the frills
```
bash scripts/run-w-frills <params-file> <profile name from nextflow.config>
```
Example
```
bash scripts/run-w-frills example_parameters.yml standard
```

# Data Processing For <a href="riboseq.org">RiboSeq.org</a>

### Automated processing of Ribo-Seq (and associated RNA-Seq) data for <a href="https://gwips.ucc.ie/">GWIPS-Viz<a> and <a href="https://trips.ucc.ie/">TRIPS-Viz</a>

## About Riboseq.org 
  
  This is a set of resources for the analysis and visualisation of publically available ribosome profiling data produced and maintained by various members of LAPTI lab in the School of Biochemistry and Cell Biology at Univeristy College Cork. These resources are well documented in their respective publications 
  
  - GWIPS-Viz
    - <a href="https://doi.org/10.1093/nar/gkx790">GWIPS-viz: 2018 update (2018).</a> Nucleic Acids Res
    - <a href="https://doi.org/10.1002/cpbi.50">The GWIPS-viz Browser (2018).</a> Current Protocols in Bioinformatics
    - <a href="http://dx.doi.org/10.1002/pmic.201400603 ">GWIPS-viz as a tool for exploring ribosome profiling evidence supporting the synthesis of alternative proteoforms (2015).</a> Proteomics
    - <a href="http://dx.doi.org/10.1093/nar/gkt1035"> GWIPS-viz: development of a ribo-seq genome browser (2014).</a> Nucleic Acids Res 
  
  - <a href="https://doi.org/10.1093/nar/gky842">Trips-Viz: a transcriptome browser for exploring Ribo-Seq data (2019).</a> Nucleic Acids Res
  
  - <a href="http://dx.doi.org/10.1080/15476286.2016.1141862">RiboGalaxy: a browser based platform for the alignment, analysis and visualization of ribosome profiling data.</a> RNA Biology-Viz
    - ** Note: Ribogalaxy is being updated currently and functionality will be restored shortly (14-2-2022)**
  #

