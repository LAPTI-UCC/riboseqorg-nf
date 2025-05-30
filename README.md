
# RiboSeqOrg Data Processing Pipeline

This Nextflow pipeline performs the processing of Ribosome Profiling (Ribo-seq) data from raw sequencing files to analysis-ready formats. The pipeline includes data acquisition, quality control, alignment, post-processing and analysis steps.


## Table of Contents
- [About RiboSeq.Org](#about-riboseqorg)
- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
  - [1. Data Acquisition](#1-data-acquisition)
  - [2. Quality Control](#2-quality-control)
  - [3. Alignment](#3-alignment)
  - [4. Post-Processing](#4-post-processing)
  - [5. Analysis](#5-analysis)
- [Parameters](#parameters)
- [Output Files](#output-files)


## About Riboseq.org 
  
  This is a set of resources for the analysis and visualisation of publically available ribosome profiling data produced and maintained by various members of LAPTI lab in the School of Biochemistry and Cell Biology at Univeristy College Cork. These resources are well documented in their respective publications 

  - RiboSeq.Org (RiboSeq Data Portal)
    - <a href="https://doi.org/10.1093/nar/gkae1020">RiboSeq.Org: an integrated suite of resources for ribosome profiling data analysis and visualization. (2025)</a> Nucleic Acids Res
  
  - GWIPS-Viz
    - <a href="https://doi.org/10.1093/nar/gkx790">GWIPS-viz: 2018 update (2018).</a> Nucleic Acids Res
    - <a href="https://doi.org/10.1002/cpbi.50">The GWIPS-viz Browser (2018).</a> Current Protocols in Bioinformatics
    - <a href="http://dx.doi.org/10.1002/pmic.201400603 ">GWIPS-viz as a tool for exploring ribosome profiling evidence supporting the synthesis of alternative proteoforms (2015).</a> Proteomics
    - <a href="http://dx.doi.org/10.1093/nar/gkt1035"> GWIPS-viz: development of a ribo-seq genome browser (2014).</a> Nucleic Acids Res 
  
  - Trips-Viz
    - <a href="https://doi.org/10.1093/nar/gky842">Trips-Viz: a transcriptome browser for exploring Ribo-Seq data (2019).</a> Nucleic Acids Res
    - <a href="https://doi.org/10.1093/nar/gkab323">Trips-Viz: an environment for the analysis of public and user-generated ribosome profiling data.</a> Nucleic Acids Res
  
  - RiboGalaxy
    - <a href="http://dx.doi.org/10.1080/15476286.2016.1141862">RiboGalaxy: a browser based platform for the alignment, analysis and visualization of ribosome profiling data.</a> RNA Biology-Viz
    - <a href="https://doi.org/10.1016/j.jmb.2023.168043">RiboGalaxy: A Galaxy-based Web Platform for Ribosome Profiling Data Processing – 2023 Update.</a> Journal of Molecular Biology


## Pipeline Overview

This pipeline processes Ribo-seq data through the following main steps:

1. **Data Acquisition**: Locates or downloads raw sequencing files
2. **Quality Control**: Validates reads for Ribosome Protected Fragment (RPF) characteristics
3. **Alignment**: Maps reads to the transcriptome and genome
4. **Post-Processing**: Converts aligned reads to usable formats (BAM, BED, BigWig, SQLite)
5. **Analysis**: Performs quality metrics and analysis with RiboMetric

## Requirements

- Nextflow (>= 21.04.0)
- Container engine: Docker, Singularity, or Conda environments
- Reference files:
  - Transcriptome and genome indices
  - Annotation SQLite file
  - Chromosome sizes file
  - GTF annotation file
  - rRNA index (for rRNA filtering)

## Usage

Basic execution of the pipeline:

```bash
nextflow run main.nf --sample_sheet [PATH_TO_SAMPLE_SHEET] -profile [PROFILE]
```

Available profiles:
- `standard`: For local execution
- `hpc_slurm`: For execution on HPC with SLURM scheduler
- `lsf`: For execution on HPC with LSF scheduler
- `docker`: Using Docker containers
- `singularity`: Using Singularity containers
- `conda`: Using Conda environments

## Pipeline Steps

### 1. Data Acquisition

**Workflow**: `DATA_ACQUISITION`

This step handles locating and retrieving the raw sequencing data.

**Processes**:
- `LOCATE`: Checks if collapsed reads already exist
  ```bash
  # Command explanation:
  # Checks if a collapsed file exists at the expected path or with _1 suffix
  collapsed_file="${params.collapsed_read_path}/${run[0..5]}/${run}.collapsed.fa.gz"
  collapsed_file_1="${params.collapsed_read_path}/${run[0..5]}/${run}_1.collapsed.fa.gz"
  
  # If file exists, create symlink; otherwise mark for processing
  if [ -f "$collapsed_file" ]; then
      ln -s "$collapsed_file" "${run}.collapsed.fa.gz"
  elif [ -f "$collapsed_file_1" ]; then
      ln -s "$collapsed_file_1" "${run}.collapsed.fa.gz"
  else
      touch "${run}_needs_processing"
  fi
  ```

- `collapse` Subworkflow: Processes raw reads if needed
  - `FASTQ_DL`: Downloads raw FASTQ files
    ```bash
    # Command:
    fastq-dl -a $run --cpus $task.cpus $args
    
    # Flag explanation:
    # -a: SRA accession number
    # --cpus: Number of CPUs to use
    # $args: Any additional arguments passed to the process
    ```
  
  - `FASTQC`: Quality assessment of raw reads
    ```bash
    # Command:
    fastqc $args --threads $task.cpus --outdir . --extract -q $fastq --adapters $adapter_list
    
    # Flag explanation:
    # --threads: Number of threads to use for processing
    # --outdir: Output directory for FastQC results
    # --extract: Extract the zipped output file
    # -q: Quiet mode
    # --adapters: Path to custom adapter list
    ```
  
  - `FIND_ADAPTERS`: Detects adapter sequences
    ```bash
    # Command:
    python3 $projectDir/bin/get_adapters.py -i $fastqc_data -a $projectDir/scripts/adapter_list.tsv -o "${prefix}_adapter_report.fa" $args
    
    # Flag explanation:
    # -i: Input FastQC data file
    # -a: Adapter list file
    # -o: Output file for detected adapters
    ```
  
  - `FASTP`: Trims adapter sequences (two steps)
    ```bash
    # Step 1: Using provided adapters
    fastp -i $reads -o ${prefix}_clipped_provided.fastq.gz $args --length_required 20 --adapter_fasta $adapter_fasta --json ${prefix}_provided_fastp.json --html ${prefix}_provided_fastp.html --thread $task.cpus
    
    # Step 2: Using auto-detection on the output from step 1
    fastp -i ${prefix}_clipped_provided.fastq.gz -o ${prefix}_clipped_final.fastq.gz $args --length_required 20 --auto_detect_adapter --json ${prefix}_auto_fastp.json --html ${prefix}_auto_fastp.html --thread $task.cpus
    
    # Flag explanation:
    # -i: Input file
    # -o: Output file
    # --length_required: Minimum read length after trimming (20 bp)
    # --adapter_fasta: FASTA file with adapter sequences
    # --auto_detect_adapter: Enable automatic adapter detection
    # --json: Output JSON report
    # --html: Output HTML report
    # --thread: Number of threads to use
    ```
  
  - `COLLAPSE_FASTQ`: Collapses duplicate reads
    ```bash
    # Command:
    RDP-Tools collapse $args $fastq
    gzip *fastq.collapsed.fa
    
    # Flag explanation:
    # collapse: RDP-Tools subcommand to collapse identical reads
    # $args: Additional arguments passed to the process
    # The output is then gzipped
    ```

**Parameters**:
- `params.fetch`: Boolean to enable/disable fetching data
- `params.force_fetch`: Boolean to force re-fetching even if files exist
- `params.collapsed_read_path`: Path to look for existing collapsed reads

### 2. Quality Control

**Workflow**: `QUALITY_CONTROL`

This step validates the reads are suitable for Ribo-seq analysis.

**Processes**:
- `GET_RPF`: Tests if reads look like genuine Ribosome Protected Fragments (RPFs)
  ```bash
  # Command:
  getRPF check ${input_file} --format collapsed --output ${prefix}_report.txt --count-pattern "${count_pattern}" $args
  
  # Flag explanation:
  # check: getRPF subcommand to perform RPF validation checks
  # --format collapsed: Specifies the input is in collapsed read format
  # --output: Path for the output report
  # --count-pattern: Regular expression pattern to extract read counts from headers
  #                 (default is typically "_x(\d+)$" to extract counts from headers like "read_x10")
  # --max-reads: (from nextflow.config) Number of reads to sample for analysis
  
  # The process validates:
  # 1. Read length distribution (typically 26-34nt for RPFs)
  # 2. 5' and 3' nucleotide preferences
  ```

**Parameters**:
- `params.collapsed_read_header_pattern`: Pattern for read headers in collapsed reads (e.g., "_x(\d+)$")
- `params.max_reads`: (from ext.args) Maximum number of reads to analyze for QC

### 3. Alignment

**Workflow**: `ALIGNMENT`

This step aligns reads to the transcriptome and genome.

**Processes**:
- `STAR_ALIGN`: Aligns reads to genome using STAR
  ```bash
  # Command:
  STAR \
      --genomeDir $index \
      --readFilesIn $reads \
      --runThreadN ${task.cpus} \
      --outFileNamePrefix ${prefix}. \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax ${params.max_multimap} \
      --outFilterMismatchNmax ${params.mismatches} \
      $output_transcriptome_bam \
      $alignment_type \
      $allow_introns \
      $trim_front \
      --readFilesCommand $unzip_command \
      --sjdbGTFfile $gtf \
      $args
      
  # Flag explanation:
  # --genomeDir: Directory containing STAR genome index
  # --readFilesIn: Input read files
  # --runThreadN: Number of threads to use
  # --outFileNamePrefix: Prefix for output files
  # --outSAMtype: Output format (BAM SortedByCoordinate)
  # --outFilterMultimapNmax: Maximum number of loci a read can map to (defaults vary)
  # --outFilterMismatchNmax: Maximum number of mismatches allowed
  # --quantMode TranscriptomeSAM: (optional) Output aligned reads in transcriptome coordinates
  # --alignEndsType Local: (optional) Performs local alignment (soft-clipping ends)
  # --alignIntronMax/--alignMatesGapMax: (optional) Max intron/gap size for alignments
  # --clip3pNbases: Number of bases to trim from the 3' end
  # --readFilesCommand: Command to decompress input files (zcat for gzipped)
  # --sjdbGTFfile: Gene annotation GTF file
  ```

- `BOWTIE_ALIGN_SORT`: Aligns reads to reference using Bowtie
  ```bash
  # Command:
  bowtie \
      -S -a --norc \
      -p ${task.cpus} \
      -v 3 \
      --seedlen 25 \
      -x $INDEX \
      -f ${reads} 2> ${meta.id}_bowtie_alignment_stats.log | samtools view -bS - | samtools sort -@ ${task.cpus} -o ${prefix}.bam
  
  samtools index ${prefix}.bam
  
  # Flag explanation:
  # -S: Output in SAM format
  # -a: Report all valid alignments per read
  # --norc: Do not align to reverse-complement reference strand
  # -p: Number of threads/cores to use
  # -v 3: Allow up to 3 mismatches in alignment
  # --seedlen 25: Use 25bp seed length for alignment
  # -x: Path to Bowtie index
  # -f: Input is in FASTA format
  # 2>: Redirect alignment statistics to log file
  
  # Then pipe to:
  # samtools view -bS: Convert SAM to BAM format
  # samtools sort: Sort the BAM file
  # -@: Number of threads for sorting
  # -o: Output file path
  
  # Finally, index the sorted BAM file
  ```

- `SAMTOOLS_INDEX`: Indexes BAM files
  ```bash
  # Command:
  samtools index -@ ${task.cpus} ${args} ${sorted_bam}
  
  # Flag explanation:
  # index: samtools subcommand to create BAM index
  # -@: Number of threads to use
  # ${args}: Any additional arguments
  # ${sorted_bam}: Path to sorted BAM file to index
  ```

**Parameters**:
- `params.star_index`: Path to STAR index
- `params.bowtie_index`: Path to Bowtie index
- `params.gtf`: Path to GTF annotation
- `params.trim_front`: Number of bases to trim from the 5' end
- `params.alignment_type`: Alignment type (Local or End-to-End)
- `params.allow_introns`: Boolean to allow introns in alignment
- `params.max_multimap`: Maximum number of multiple mappings to output
- `params.mismatches`: Maximum number of mismatches allowed
- `params.save_star_transcriptome_bam`: Boolean to generate transcriptome BAM

### 4. Post-Processing

**Workflow**: `POST_PROCESSING`

This step processes aligned reads into different formats for analysis and visualization.

**Processes**:
- `SAMTOOLS_NAME_SORT`: Sorts BAM by read name
  ```bash
  # Command:
  samtools sort \
      -n \
      -@ ${task.cpus} \
      ${args} \
      -o ${prefix}.name_sorted.bam \
      ${transcriptome_alignments}
  
  # Flag explanation:
  # sort: samtools subcommand for sorting
  # -n: Sort by read name (instead of coordinate)
  # -@: Number of threads to use
  # ${args}: Additional arguments
  # -o: Output file path
  # ${transcriptome_alignments}: Input BAM file
  
  # Note: Name sorting is required for bam_to_sqlite.py, which processes reads by name
  ```

- `BAM_TO_SQLITE`: Converts BAM to SQLite database
  ```bash
  # Command:
  python3 $projectDir/bin/bam_to_sqlite.py \
      --bam ${sorted_bam} \
      --annotation ${annotation_sqlite} \
      --output ${prefix}.sqlite \
      $args
  
  # Flag explanation:
  # --bam: Input BAM file (must be name-sorted)
  # --annotation: SQLite annotation database containing transcript information
  # --output: Output SQLite file path
  
  # This process:
  # 1. Reads alignments from BAM file
  # 2. Maps reads to transcript positions
  # 3. Calculates A-site positions using read lengths
  # 4. Stores read counts by position, read length, and ambiguity
  # 5. Calculates triplet periodicity and other metrics
  # 6. Determines optimal A-site offsets for each read length
  ```

- `EXTRACT_OFFSETS`: Extracts A-site offsets from SQLite
  ```bash
  # Command:
  python $projectDir/bin/extract_offsets.py \
      -d ${sqlite_file} \
      -o ${prefix}.offsets.txt \
      $args
  
  # Flag explanation:
  # -d: Input SQLite file
  # -o: Output text file for offsets
  # $args: May include '-e fiveprime' or '-e threeprime' to specify which end to extract offsets for
  
  # This script extracts the optimal A-site offsets for each read length
  # Output format is a tab-delimited file with columns:
  # read_length    offset
  ```

- `BAM_TO_BED`: Converts BAM to BED format with A-site mapping
  ```bash
  # Command:
  python3 $projectDir/bin/bam_to_bed.py \
      --bam ${bam} \
      --offsets ${offsets} \
      --prefix ${prefix} \
      $args
  
  sort -k1,1 -k2,2n ${prefix}.bedgraph > ${prefix}.sorted.bedgraph
  
  # Flag explanation:
  # --bam: Input BAM file
  # --offsets: Text file with read length-specific A-site offsets
  # --prefix: Output prefix
  # $args: Additional arguments, possibly including '--stranded' to generate strand-specific output
  
  # The script:
  # 1. Reads offset file to get length-specific A-site offsets
  # 2. Maps each read to the A-site position based on read length
  # 3. Outputs coverage in BedGraph format
  # 4. Can generate strand-specific outputs if --stranded is specified
  
  # The output is then sorted using UNIX sort:
  # -k1,1: Sort by first column (chromosome)
  # -k2,2n: Then sort numerically by second column (position)
  ```

- `BEDGRAPH_TO_BIGWIG`: Converts BedGraph to BigWig
  ```bash
  # Command:
  bedGraphToBigWig ${bedgraph} ${chrom_sizes} ${prefix}.bw
  
  # Parameters:
  # ${bedgraph}: Input sorted BedGraph file
  # ${chrom_sizes}: File with chromosome sizes
  # ${prefix}.bw: Output BigWig file
  
  # This tool converts the text-based BedGraph format to the binary BigWig format
  # BigWig files are more efficient for visualization in genome browsers
  ```

**Parameters**:
- `params.annotation_sqlite`: Path to SQLite annotation database
- `params.chrom_sizes_file`: Path to chromosome sizes file

### 5. Analysis

**Workflow**: `ANALYSIS`

This step performs quality assessment and analysis of the processed data.

**Processes**:
- `RIBOMETRIC`: Runs RiboMetric on the processed data
  ```bash
  # Command:
  RiboMetric run \
      --bam ${transcriptome_bam} \
      --annotation ${ribometric_annotation} \
      --threads $task.cpus \
      --html \
      --json \
      --csv \
      --offset-calculation-method tripsviz \
      -S 10000000
      $args
  
  # Flag explanation:
  # run: RiboMetric subcommand to perform analysis
  # --bam: Input BAM file with transcriptome alignments
  # --annotation: Path to RiboMetric annotation file
  # --threads: Number of CPU threads to use
  # --html: Generate HTML report
  # --json: Generate JSON output
  # --csv: Generate CSV output
  # --offset-calculation-method tripsviz: Use TripsViz method to calculate A-site offsets
  # -S 10000000: Subsample to 10 million reads (for performance)
  
  # RiboMetric analyzes:
  # 1. Read length distribution
  # 2. Metagene profiles around start and stop codons
  # 3. 3-nucleotide periodicity
  # 4. P-site/A-site offset calculation
  # 5. Codon usage
  # 6. Reading frame distribution
  # 7. Coverage characteristics
  ```

**Parameters**:
- `params.ribometric_annotation`: Path to RiboMetric annotation (SQLite format)

**Note**: The RiboMetric analysis is optional in the main workflow in main.nf, as its invocation is commented out in the most recent version.

## Parameters

The pipeline accepts the following main parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `params.sample_sheet` | Path to sample sheet CSV | - |
| `params.outdir` | Output directory | - |
| `params.fetch` | Enable data download | `true` |
| `params.force_fetch` | Force re-download | `false` |
| `params.star_index` | Path to STAR index | - |
| `params.bowtie_index` | Path to Bowtie index | - |
| `params.gtf` | Path to GTF file | - |
| `params.annotation_sqlite` | Path to annotation SQLite | - |
| `params.chrom_sizes_file` | Path to chromosome sizes | - |
| `params.rRNA_index` | Path to rRNA index | - |
| `params.collapsed_read_header_pattern` | Pattern for read headers | `_x(\d+)# RiboSeqOrg Data Processing Pipeline

This Nextflow pipeline performs the processing of Ribosome Profiling (Ribo-seq) data from raw sequencing files to analysis-ready formats. The pipeline includes data acquisition, quality control, alignment, post-processing and analysis steps.

## Table of Contents
- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
  - [1. Data Acquisition](#1-data-acquisition)
  - [2. Quality Control](#2-quality-control)
  - [3. Alignment](#3-alignment)
  - [4. Post-Processing](#4-post-processing)
  - [5. Analysis](#5-analysis)
- [Parameters](#parameters)
- [Output Files](#output-files)

## Pipeline Overview

This pipeline processes Ribo-seq data through the following main steps:

1. **Data Acquisition**: Locates or downloads raw sequencing files
2. **Quality Control**: Validates reads for Ribosome Protected Fragment (RPF) characteristics
3. **Alignment**: Maps reads to the transcriptome and genome
4. **Post-Processing**: Converts aligned reads to usable formats (BAM, BED, BigWig, SQLite)
5. **Analysis**: Performs quality metrics and analysis with RiboMetric

## Requirements

- Nextflow (>= 21.04.0)
- Container engine: Docker, Singularity, or Conda environments
- Reference files:
  - Transcriptome and genome indices
  - Annotation SQLite file
  - Chromosome sizes file
  - GTF annotation file
  - rRNA index (for rRNA filtering)

## Usage

Basic execution of the pipeline:

```bash
nextflow run main.nf --sample_sheet [PATH_TO_SAMPLE_SHEET] -profile [PROFILE]
```

Available profiles:
- `standard`: For local execution
- `hpc_slurm`: For execution on HPC with SLURM scheduler
- `lsf`: For execution on HPC with LSF scheduler
- `docker`: Using Docker containers
- `singularity`: Using Singularity containers
- `conda`: Using Conda environments

## Pipeline Steps

### 1. Data Acquisition

**Workflow**: `DATA_ACQUISITION`

This step handles locating and retrieving the raw sequencing data.

**Processes**:
- `LOCATE`: Checks if collapsed reads already exist
  ```bash
  # Command explanation:
  # Checks if a collapsed file exists at the expected path or with _1 suffix
  collapsed_file="${params.collapsed_read_path}/${run[0..5]}/${run}.collapsed.fa.gz"
  collapsed_file_1="${params.collapsed_read_path}/${run[0..5]}/${run}_1.collapsed.fa.gz"
  
  # If file exists, create symlink; otherwise mark for processing
  if [ -f "$collapsed_file" ]; then
      ln -s "$collapsed_file" "${run}.collapsed.fa.gz"
  elif [ -f "$collapsed_file_1" ]; then
      ln -s "$collapsed_file_1" "${run}.collapsed.fa.gz"
  else
      touch "${run}_needs_processing"
  fi
  ```

- `collapse` Subworkflow: Processes raw reads if needed
  - `FASTQ_DL`: Downloads raw FASTQ files
    ```bash
    # Command:
    fastq-dl -a $run --cpus $task.cpus $args
    
    # Flag explanation:
    # -a: SRA accession number
    # --cpus: Number of CPUs to use
    # $args: Any additional arguments passed to the process
    ```
  
  - `FASTQC`: Quality assessment of raw reads
    ```bash
    # Command:
    fastqc $args --threads $task.cpus --outdir . --extract -q $fastq --adapters $adapter_list
    
    # Flag explanation:
    # --threads: Number of threads to use for processing
    # --outdir: Output directory for FastQC results
    # --extract: Extract the zipped output file
    # -q: Quiet mode
    # --adapters: Path to custom adapter list
    ```
  
  - `FIND_ADAPTERS`: Detects adapter sequences
    ```bash
    # Command:
    python3 $projectDir/bin/get_adapters.py -i $fastqc_data -a $projectDir/scripts/adapter_list.tsv -o "${prefix}_adapter_report.fa" $args
    
    # Flag explanation:
    # -i: Input FastQC data file
    # -a: Adapter list file
    # -o: Output file for detected adapters
    ```
  
  - `FASTP`: Trims adapter sequences (two steps)
    ```bash
    # Step 1: Using provided adapters
    fastp -i $reads -o ${prefix}_clipped_provided.fastq.gz $args --length_required 20 --adapter_fasta $adapter_fasta --json ${prefix}_provided_fastp.json --html ${prefix}_provided_fastp.html --thread $task.cpus
    
    # Step 2: Using auto-detection on the output from step 1
    fastp -i ${prefix}_clipped_provided.fastq.gz -o ${prefix}_clipped_final.fastq.gz $args --length_required 20 --auto_detect_adapter --json ${prefix}_auto_fastp.json --html ${prefix}_auto_fastp.html --thread $task.cpus
    
    # Flag explanation:
    # -i: Input file
    # -o: Output file
    # --length_required: Minimum read length after trimming (20 bp)
    # --adapter_fasta: FASTA file with adapter sequences
    # --auto_detect_adapter: Enable automatic adapter detection
    # --json: Output JSON report
    # --html: Output HTML report
    # --thread: Number of threads to use
    ```
  
  - `COLLAPSE_FASTQ`: Collapses duplicate reads
    ```bash
    # Command:
    RDP-Tools collapse $args $fastq
    gzip *fastq.collapsed.fa
    
    # Flag explanation:
    # collapse: RDP-Tools subcommand to collapse identical reads
    # $args: Additional arguments passed to the process
    # The output is then gzipped
    ```

**Parameters**:
- `params.fetch`: Boolean to enable/disable fetching data
- `params.force_fetch`: Boolean to force re-fetching even if files exist
- `params.collapsed_read_path`: Path to look for existing collapsed reads

### 2. Quality Control

**Workflow**: `QUALITY_CONTROL`

This step validates the reads are suitable for Ribo-seq analysis.

**Processes**:
- `GET_RPF`: Tests if reads look like genuine Ribosome Protected Fragments (RPFs)
  ```bash
  # Command:
  getRPF check ${input_file} --format collapsed --output ${prefix}_report.txt --count-pattern "${count_pattern}" --max-reads 100000

  # Flag explanation:
  # check: getRPF subcommand to perform RPF validation checks
  # --format collapsed: Specifies the input is in collapsed read format
  # --output: Path for the output report
  # --count-pattern: Regular expression pattern to extract read counts from headers (typically "_x(\d+)$")
  # --max-reads: Maximum number of reads to sample for analysis (default: 100000 in your pipeline)
  ```
#### RPF-Specific Checks Performed:

**Length Distribution Check**:
RPFs typically have a narrow length distribution between 26-35 nucleotides
- `PASS`: Majority of reads (>50%) fall within expected length range, with peak length in range
- `WARN`: Low fraction (30-50%) of reads in RPF range, or peak length outside range
- `FAIL`: Very low fraction (<30%) of reads in expected RPF range


**Base Composition Check**:
Examines positions with extreme nucleotide bias
- `PASS`: No positions have a single nucleotide frequency >85%
- `FAIL`: One or more positions have extreme bias (>85% one nucleotide)
Identifies adapter contamination or sequence artifacts


**GC Content Check**:
- `PASS`: GC content between 35-65%
- `WARN`: GC content outside expected range
Abnormal GC can indicate contamination or bias

**Parameters**:
- `params.collapsed_read_header_pattern`: Pattern for read headers in collapsed reads (e.g., "_x(\d+)$")
- `params.max_reads`: (from ext.args) Maximum number of reads to analyze for QC

### 3. Alignment

**Workflow**: `ALIGNMENT`

This step aligns reads to the transcriptome and genome.

**Processes**:
- `STAR_ALIGN`: Aligns reads to genome using STAR
  ```bash
  # Command:
  STAR \
      --genomeDir $index \
      --readFilesIn $reads \
      --runThreadN ${task.cpus} \
      --outFileNamePrefix ${prefix}. \
      --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax ${params.max_multimap} \
      --outFilterMismatchNmax ${params.mismatches} \
      $output_transcriptome_bam \
      $alignment_type \
      $allow_introns \
      $trim_front \
      --readFilesCommand $unzip_command \
      --sjdbGTFfile $gtf \
      $args
      
  # Flag explanation:
  # --genomeDir: Directory containing STAR genome index
  # --readFilesIn: Input read files
  # --runThreadN: Number of threads to use
  # --outFileNamePrefix: Prefix for output files
  # --outSAMtype: Output format (BAM SortedByCoordinate)
  # --outFilterMultimapNmax: Maximum number of loci a read can map to (defaults vary)
  # --outFilterMismatchNmax: Maximum number of mismatches allowed
  # --quantMode TranscriptomeSAM: (optional) Output aligned reads in transcriptome coordinates
  # --alignEndsType Local: (optional) Performs local alignment (soft-clipping ends)
  # --alignIntronMax/--alignMatesGapMax: (optional) Max intron/gap size for alignments
  # --clip3pNbases: Number of bases to trim from the 3' end
  # --readFilesCommand: Command to decompress input files (zcat for gzipped)
  # --sjdbGTFfile: Gene annotation GTF file
  ```

- `BOWTIE_ALIGN_SORT`: Aligns reads to reference using Bowtie
  ```bash
  # Command:
  bowtie \
      -S -a --norc \
      -p ${task.cpus} \
      -v 3 \
      --seedlen 25 \
      -x $INDEX \
      -f ${reads} 2> ${meta.id}_bowtie_alignment_stats.log | samtools view -bS - | samtools sort -@ ${task.cpus} -o ${prefix}.bam
  
  samtools index ${prefix}.bam
  
  # Flag explanation:
  # -S: Output in SAM format
  # -a: Report all valid alignments per read
  # --norc: Do not align to reverse-complement reference strand
  # -p: Number of threads/cores to use
  # -v 3: Allow up to 3 mismatches in alignment
  # --seedlen 25: Use 25bp seed length for alignment
  # -x: Path to Bowtie index
  # -f: Input is in FASTA format
  # 2>: Redirect alignment statistics to log file
  
  # Then pipe to:
  # samtools view -bS: Convert SAM to BAM format
  # samtools sort: Sort the BAM file
  # -@: Number of threads for sorting
  # -o: Output file path
  
  # Finally, index the sorted BAM file
  ```

- `SAMTOOLS_INDEX`: Indexes BAM files
  ```bash
  # Command:
  samtools index -@ ${task.cpus} ${args} ${sorted_bam}
  
  # Flag explanation:
  # index: samtools subcommand to create BAM index
  # -@: Number of threads to use
  # ${args}: Any additional arguments
  # ${sorted_bam}: Path to sorted BAM file to index
  ```

**Parameters**:
- `params.star_index`: Path to STAR index
- `params.bowtie_index`: Path to Bowtie index
- `params.gtf`: Path to GTF annotation
- `params.trim_front`: Number of bases to trim from the 5' end
- `params.alignment_type`: Alignment type (Local or End-to-End)
- `params.allow_introns`: Boolean to allow introns in alignment
- `params.max_multimap`: Maximum number of multiple mappings to output
- `params.mismatches`: Maximum number of mismatches allowed
- `params.save_star_transcriptome_bam`: Boolean to generate transcriptome BAM

### 4. Post-Processing

**Workflow**: `POST_PROCESSING`

This step processes aligned reads into different formats for analysis and visualization.

**Processes**:
- `SAMTOOLS_NAME_SORT`: Sorts BAM by read name
  ```bash
  # Command:
  samtools sort \
      -n \
      -@ ${task.cpus} \
      ${args} \
      -o ${prefix}.name_sorted.bam \
      ${transcriptome_alignments}
  
  # Flag explanation:
  # sort: samtools subcommand for sorting
  # -n: Sort by read name (instead of coordinate)
  # -@: Number of threads to use
  # ${args}: Additional arguments
  # -o: Output file path
  # ${transcriptome_alignments}: Input BAM file
  
  # Note: Name sorting is required for bam_to_sqlite.py, which processes reads by name
  ```

- `BAM_TO_SQLITE`: Converts BAM to SQLite database
  ```bash
  # Command:
  python3 $projectDir/bin/bam_to_sqlite.py \
      --bam ${sorted_bam} \
      --annotation ${annotation_sqlite} \
      --output ${prefix}.sqlite \
      $args
  
  # Flag explanation:
  # --bam: Input BAM file (must be name-sorted)
  # --annotation: SQLite annotation database containing transcript information
  # --output: Output SQLite file path
  
  # This process:
  # 1. Reads alignments from BAM file
  # 2. Maps reads to transcript positions
  # 3. Calculates A-site positions using read lengths
  # 4. Stores read counts by position, read length, and ambiguity
  # 5. Calculates triplet periodicity and other metrics
  # 6. Determines optimal A-site offsets for each read length
  ```

- `EXTRACT_OFFSETS`: Extracts A-site offsets from SQLite
  ```bash
  # Command:
  python $projectDir/bin/extract_offsets.py \
      -d ${sqlite_file} \
      -o ${prefix}.offsets.txt \
      $args
  
  # Flag explanation:
  # -d: Input SQLite file
  # -o: Output text file for offsets
  # $args: May include '-e fiveprime' or '-e threeprime' to specify which end to extract offsets for
  
  # This script extracts the optimal A-site offsets for each read length
  # Output format is a tab-delimited file with columns:
  # read_length    offset
  ```

- `BAM_TO_BED`: Converts BAM to BED format with A-site mapping
  ```bash
  # Command:
  python3 $projectDir/bin/bam_to_bed.py \
      --bam ${bam} \
      --offsets ${offsets} \
      --prefix ${prefix} \
      $args
  
  sort -k1,1 -k2,2n ${prefix}.bedgraph > ${prefix}.sorted.bedgraph
  
  # Flag explanation:
  # --bam: Input BAM file
  # --offsets: Text file with read length-specific A-site offsets
  # --prefix: Output prefix
  # $args: Additional arguments, possibly including '--stranded' to generate strand-specific output
  
  # The script:
  # 1. Reads offset file to get length-specific A-site offsets
  # 2. Maps each read to the A-site position based on read length
  # 3. Outputs coverage in BedGraph format
  # 4. Can generate strand-specific outputs if --stranded is specified
  
  # The output is then sorted using UNIX sort:
  # -k1,1: Sort by first column (chromosome)
  # -k2,2n: Then sort numerically by second column (position)
  ```

- `BEDGRAPH_TO_BIGWIG`: Converts BedGraph to BigWig
  ```bash
  # Command:
  bedGraphToBigWig ${bedgraph} ${chrom_sizes} ${prefix}.bw
  
  # Parameters:
  # ${bedgraph}: Input sorted BedGraph file
  # ${chrom_sizes}: File with chromosome sizes
  # ${prefix}.bw: Output BigWig file
  
  # This tool converts the text-based BedGraph format to the binary BigWig format
  # BigWig files are more efficient for visualization in genome browsers
  ```

**Parameters**:
- `params.annotation_sqlite`: Path to SQLite annotation database
- `params.chrom_sizes_file`: Path to chromosome sizes file

### 5. Analysis

**Workflow**: `ANALYSIS`

This step performs quality assessment and analysis of the processed data.

**Processes**:
- `RIBOMETRIC`: Runs RiboMetric on the processed data
  ```bash
  # Command:
  RiboMetric run \
      --bam ${transcriptome_bam} \
      --annotation ${ribometric_annotation} \
      --threads $task.cpus \
      --html \
      --json \
      --csv \
      --offset-calculation-method tripsviz \
      -S 10000000
      $args
  
  # Flag explanation:
  # run: RiboMetric subcommand to perform analysis
  # --bam: Input BAM file with transcriptome alignments
  # --annotation: Path to RiboMetric annotation file
  # --threads: Number of CPU threads to use
  # --html: Generate HTML report
  # --json: Generate JSON output
  # --csv: Generate CSV output
  # --offset-calculation-method tripsviz: Use TripsViz method to calculate A-site offsets
  # -S 10000000: Subsample to 10 million reads (for performance)
  
  # RiboMetric analyzes:
  # 1. Read length distribution
  # 2. Metagene profiles around start and stop codons
  # 3. 3-nucleotide periodicity
  # 4. P-site/A-site offset calculation
  # 5. Codon usage
  # 6. Reading frame distribution
  # 7. Coverage characteristics
  ```

**Parameters**:
- `params.ribometric_annotation`: Path to RiboMetric annotation (SQLite format)

**Note**: The RiboMetric analysis is optional in the main workflow in main.nf, as its invocation is commented out in the most recent version.

 |
| `params.collapsed_read_path` | Path to look for existing collapsed reads | - |
| `params.trim_front` | Number of bases to trim from 5' end | `0` |
| `params.alignment_type` | Alignment type (Local or End-to-End) | - |
| `params.allow_introns` | Allow introns in alignment | - |
| `params.max_multimap` | Maximum number of multiple mappings | - |
| `params.mismatches` | Maximum number of mismatches allowed | - |
| `params.save_star_transcriptome_bam` | Generate transcriptome BAM from STAR | `false` |
| `params.ribometric_annotation` | Path to RiboMetric annotation | - |

### Sample Sheet Format

The pipeline expects a sample sheet CSV with the following columns:

```
SRAStudy,Run,ScientificName
```

Example:
```
SRAStudy,Run,ScientificName
SRP123456,SRR1234567,Homo sapiens
SRP123456,SRR1234568,Homo sapiens
```

### Configuration Profiles

The pipeline includes several configuration profiles, specified with `-profile`:

| Profile | Description |
|---------|-------------|
| `standard` | Default configuration for local execution |
| `hpc_slurm` | Configuration for HPC clusters using SLURM scheduler |
| `lsf` | Configuration for HPC clusters using LSF scheduler |
| `docker` | Uses Docker containers for tools |
| `singularity` | Uses Singularity containers for tools |
| `conda` | Uses Conda environments for tools |

### Resource Configuration

Resource allocation for processes can be found in `nextflow.config`:

| Process Label | Memory | CPUs |
|---------------|--------|------|
| `high` | 8GB | 8 |
| `low` | 2GB | 2 |
| Others | Variable | Variable |

Some processes have specific settings:
- `FASTQ_DL`: Limited to 3 concurrent processes
- `STAR_ALIGN`: Limited to 5 concurrent processes

## Output Files

The pipeline produces the following main output directories:

| Directory | Description | File Types |
|-----------|-------------|------------|
| `fastq` | Downloaded FASTQ files | `.fastq.gz` |
| `collapsed_fa` | Collapsed FASTA files | `.fa.gz` |
| `fastqc` | FastQC quality reports | `.html`, `.zip` |
| `adapter_reports` | Detected adapter sequences | `.fa` |
| `fastp` | Trimming reports and trimmed reads | `.json`, `.html`, `_clipped_final.fastq.gz` |
| `star_align/bam` | STAR genome alignments | `.bam`, `.bai` |
| `star_align/transcriptome_bam` | STAR transcriptome alignments | `.toTranscriptome.out.bam` |
| `star_align/logs` | STAR alignment logs | `.Log.final.out` |
| `offsets` | A-site offsets from ribosome profiling | `.offsets.txt` |
| `bedgraphs` | BedGraph format coverage files | `.bedgraph`, `.sorted.bedgraph` |
| `bigwigs` | BigWig format coverage files | `.bw` |
| `sqlites` | SQLite databases with ribosome profiling data | `.sqlite` |
| `RiboMetric` | RiboMetric analysis results | `.html`, `.json`, `.csv` |
| `multiqc` | MultiQC reports aggregating QC info | `multiqc_report.html`, `_data/` |

### Key Output Files Detail

#### Offsets Files (`.offsets.txt`)
Tab-delimited text files containing read length-specific A-site offsets used for ribosome profiling analysis:
```
read_length    offset
25             12
26             13
27             13
28             14
...
```

#### SQLite Databases (`.sqlite`)
Contain extensive information about ribosome profiling data including:
- Read counts by transcript position and read length
- A-site offsets for each read length
- Translation periodicity metrics
- Read nucleotide composition
- Metagene profiles around start and stop codons

#### BedGraph Files (`.bedgraph`)
Text files with genome coverage at base resolution:
```
chromosome  start  end  coverage
chr1        1000   1001 15
chr1        1001   1002 18
...
```

#### BigWig Files (`.bw`)
Binary compressed format of the BedGraph data for efficient genome browser visualization.

