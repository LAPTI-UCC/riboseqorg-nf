# Data Processing For <a href="riboseq.org">RiboSeq.org<a>

### Automated processing of Ribo-Seq (and associated RNA-Seq) data for <a href="https://gwips.ucc.ie/">GWIPS-Viz<a> and <a href="https://trips.ucc.ie/">TRIPS-Viz<a>

## About Riboseq.org 
  
  This is a set of resources for the analysis and visualisation of publically available ribosome profiling data produced and maintained by various members of LAPTI lab in the School of Biochemistry and Cell Biology at Univeristy College Cork. These resources are well documented in their respective publications 
  
  - GWIPS-Viz
    - <a href="https://doi.org/10.1093/nar/gkx790">GWIPS-viz: 2018 update (2018).<a> Nucleic Acids Res
    - <a href="https://doi.org/10.1002/cpbi.50">The GWIPS-viz Browser (2018).<a> Current Protocols in Bioinformatics
    - <a href="http://dx.doi.org/10.1002/pmic.201400603 ">GWIPS-viz as a tool for exploring ribosome profiling evidence supporting the synthesis of alternative proteoforms (2015).<a> Proteomics
    - <a href="http://dx.doi.org/10.1093/nar/gkt1035"> GWIPS-viz: development of a ribo-seq genome browser (2014).<a> Nucleic Acids Res 
  
  - <a href="https://doi.org/10.1093/nar/gky842">Trips-Viz: a transcriptome browser for exploring Ribo-Seq data (2019).<a> Nucleic Acids Res
  
  - <a href="http://dx.doi.org/10.1080/15476286.2016.1141862">RiboGalaxy: a browser based platform for the alignment, analysis and visualization of ribosome profiling data.<a> RNA Biology-Viz
    - ** Note: Ribogalaxy is being updated currently and functionality will be restored shortly (14-2-2022)**
  #

## Requirements 
  
  - <a href="https://pypi.org/project/biopython/" > Biopython <a>
  - <a href="https://pypi.org/project/pandas/"> Pandas <a>
  - <a href="https://pypi.org/project/validators/"> Validators <a>
  
  
  #
  
  ## Outline 
  1. Produce Database Of All Available Ribosome Profiling Studies 
  2. Gather Metadata 
  3. Fetch Files and Infer Gaps in Metadata
  4. Run Pipeline 
  5. Upload to GWIPS & TRIPS
  
  # 
  
  ##  1. Produce Database Of All Available Ribosome Profiling Studies
  
  In recent years the rate at which ribosome profiling studies have been published has steadily increased. When the riboseq.org resources were initiatlly developed the number of available ribo-seq datasets was managable via manual inclusion. Here we put in place a method that records the details of relevant ribosome profiling data deposited in GEO 
  
  Initially manual searching of GEO and SRA were used along with <a href="10.3390/biology10101026">ARGEOS<a>. The outputs of each of these methods were colated to find the set of unique datasets.  
  
  
## 2. Gather Metadata
  
  GEO and SRA run tables contain valuable metadata that may be important for the processing and cateloging of the datasets. In this step we use python scripts to glean what we can from the information available 
  
## 3. Fetch Files and Infer Gaps in Metadata
  A common problem with reprocessing data for these resources is that the data is deposited in GEO and SRA with inconsistent metadata. In the stage of the process we carry out a number of steps to check for the relevant data in the provided metadata and where it is absent we infer it from the data itself. This relates to information such as cell type and treatment but also UMI position and adapter position/sequence. 
  
## 4. Run pipeline
  
  In this stage we use nextflow to process the fetched reads following the schema below
  
  ![Deptiction of the data processing pipeline](https://github.com/JackCurragh/riboseq_data_processing/blob/main/images/pipeline.drawio.png)
## 5. Upload to GWIPS and TRIPS
  
  This stage uses the metadata to upload the processed files to the web resources in an automated fashion
