#python /home/DATA/GWIPS_viz/python_and_other_scripts/Converting_to_ribosomeprofiles/bowtieOutputToRibosomeProfileBedfile_RNase1.py /home/DATA/GWIPS_viz/Ribo_seq/Gao14_human/QTI_HEK293_AminoAcid_Starvation/QTI_HEK293_AminoAcid_Starvation.bam_sorted.bam 12 /home/DATA/GWIPS_viz/Annotations_genomes_etc/Homo_sapiens_UCSC/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa


from Bio import  SeqIO
import pysam, os
import argparse

def run_weight_centered(all_reads):
	'''
	calulate asite positions using weight centered approach
	'''
	sequence = {}
		
	for read in all_reads :

		if read.qlen < 25 : continue 

		protect_nts = sorted(read.positions)

		for Asite in protect_nts[12:-12] :
			if Asite in sequence:
				sequence[Asite] += 1.0/len(protect_nts[12:-12])
			else:
				sequence[Asite] = 1.0/len(protect_nts[12:-12])

		return sequence


def run_offset(all_reads, offset):
	'''
	calculate a-site position using provided offset. 
	'''
	sequence = {}
		
	for read in all_reads :
		if read.qlen < 25 : continue

		protect_nts = sorted(read.positions)
		if "_x" in read.qname:
			read_count = int(read.qname.split("_x")[1])
			print(read.qname, read_count)
		else:
			read_count = 1

		if not read.is_reverse:
			Asite = protect_nts[offset]
		else :
			Asite = protect_nts[-1 - offset]

		if Asite in sequence:
			sequence[Asite] += read_count
		else:
			sequence[Asite] = read_count
	
	return sequence


def create_chrom_sizes(fasta_path):
	'''
	create chrom.sizes file from fasta input
	'''

	chromSizesoutput = open(fasta_path + "_chrom.sizes","w")

	records = []
	record = False
	for line in open(fasta_path, 'r').readlines():
		if line[0] == '>':
			if record:
				records.append(record)
			record = [line.strip("\n").split(' ')[0][1:], 0]

		else:
			sequence = line.strip('\n')
			record[1] += len(sequence)
			
	for seq_record in records:
		output_line = '%s\t%i\n' % (seq_record[0], seq_record[1])
		chromSizesoutput.write(output_line)

	chromSizesoutput.close()


def main(bam_path, fasta_path, mode="offset", offset=0, create_bw=False):
	'''
	create sorted bed file of a ribosome profile in two mode options
	offset - calculate a site with an inputted estimated distance from read end. Same applied to all read lengths 
	weight - weight centered approach is used (suitable for mnase digested reads)
	'''
	alignments	= pysam.Samfile(bam_path, 'rb') 

	with open(fasta_path) as in_seq_handle:
		seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
		seq_dict_keys =  sorted(seq_dict.keys())

	bed_path = str(bam_path) + ".bed"

	bedfile = open(bed_path, "w")

	for chrom in seq_dict_keys:		
		try:
			all_reads = alignments.fetch(chrom)
		except:
			print(f"No reads fetchable from provided bam file for chromosome {chrom}. Perhaps the fasta file is not what you aligned to")

		if mode == "offset":
			sequence = run_offset(all_reads, offset)
		
		elif mode == "weight":
			sequence = run_weight_centered(all_reads)

		for Asite in sorted(sequence):
			bedfile.write(f"{chrom}\t{Asite}\t{Asite + 1}\t{sequence[Asite]}\n")
			
	bedfile.close()

	command = "sort -k1,1 -k2,2n %s > %s"%(bed_path, bed_path + ".sorted")
	os.system(command)


if __name__ == '__main__':
    # Use the argparse module to parse the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_path', help='path to the BAM file')
    parser.add_argument('offset', type=int, help='offset for aligning the ribosome footprints (applied to all read lengths)')
    parser.add_argument('fasta_path', help='path to the FASTA file')
    parser.add_argument('mode', help='\'offset\' or \'weight\'')
    args = parser.parse_args()

    # Check if the BAM file is indexed, and create the index if necessary
    if not os.path.exists(args.bam_path + ".bai"):
        print(args.bam_path)
        pysam.index(args.bam_path)

    # Call the main function to process the BAM LoCAM and create the BigWig file
    main(bam_path=args.bam_path, fasta_path=args.fasta_path, mode=args.mode, offset=args.offset)



