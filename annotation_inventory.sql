CREATE TABLE annotation_inventory(
   organism varchar(100) PRIMARY KEY,
   rRNA_index varchar(100),
   transcriptome_index varchar(100),
   genome_index varchar(100),
   genome_fasta varchar(100),
   annotation_sqlite varchar(100),
   chrom_sizes_file varchar(100),
   scientific_name varchar(100)
);

CREATE TABLE primary_organism(
    scientific_name varchar(100)PRIMARY KEY,
    organism varchar(100)
);