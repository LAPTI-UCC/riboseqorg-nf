/* THE pipeline */


/* -------------------
PRE-PROCESSING BRANCH
--------------------- */



project_dir = projectDir  /// specify a new variable, the project directory ///

include { CLIP_FASTQ; rRNA_MAPPING; FASTQC_ON_PROCESSED; MULTIQC_ON_FASTQ } from "./modules/processing-tasks.nf"
include { TRANSCRIPTOME_MAPPING; TRANSCRIPTOME_SAM_TO_BAM; BAM_TO_SQLITE; RIBO_QC } from "./modules/processing-tasks.nf"
include { GENOME_MAPPING; INDEX_BAM; BAM_TO_COVBED; GENOME_BAM_TO_BED; BED_TO_BIGWIG} from "./modules/processing-tasks.nf"
include { BED_TO_BIGWIG as BED_TO_COV_BIGWIG } from "./modules/processing-tasks.nf"



workflow {

	fastq_data = Channel.fromPath ( params.fastq_files )
	CLIP_FASTQ          ( fastq_data )
	rRNA_MAPPING        ( CLIP_FASTQ.out )
	FASTQC_ON_PROCESSED ( rRNA_MAPPING.out.fastq_less_rRNA )
	MULTIQC_ON_FASTQ    ( FASTQC_ON_PROCESSED.out.fastqc_full_reports )		

    /// TRANSCRIPTOME MAPPING ///
	if ( params.skip_trips == false ) {
		
		TRANSCRIPTOME_MAPPING    ( rRNA_MAPPING.out.fastq_less_rRNA )
		TRANSCRIPTOME_SAM_TO_BAM ( TRANSCRIPTOME_MAPPING.out.transcriptome_sam )
		BAM_TO_SQLITE            ( TRANSCRIPTOME_SAM_TO_BAM.out )
		RIBO_QC					 ( BAM_TO_SQLITE.out, FASTQC_ON_PROCESSED.out.fastqc_data )

	}

    /// GENOME MAPPING ///
	if ( params.skip_gwips == false ) {

		GENOME_MAPPING        ( rRNA_MAPPING.out.fastq_less_rRNA )
		/// This block is for RNA-Seq studies only. It's executed depending on a parameter, which defines the type of study we are working with.
		params.x = "Is_RNA-Seq"
		if (params.x != "Is_RNA-Seq") {
			INDEX_BAM   	  ( GENOME_MAPPING.out.genome_sorted_bam )
			GENOME_BAM_TO_BED ( INDEX_BAM.out.genome_index_sorted_bam, INDEX_BAM.out.genome_index_sorted_bam_bai )
			BED_TO_BIGWIG     ( GENOME_BAM_TO_BED.out.sorted_beds )
		}
		BAM_TO_COVBED     ( GENOME_MAPPING.out.genome_sorted_bam )
		BED_TO_COV_BIGWIG ( BAM_TO_COVBED.out.coverage_beds )
		
    }
}