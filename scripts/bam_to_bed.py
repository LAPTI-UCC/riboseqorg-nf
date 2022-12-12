#!/usr/bin/env python3

import argparse
import pysam 
from Bio import SeqIO


def get_fasta_as_dictionary(genomic_fasta):
	'''
	Reads a genomic fasta file as a dictionary. Chromosomes are keys, sequences are values.
	'''

	in_seq_handle = open(genomic_fasta)
	seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
	in_seq_handle.close()

	return seq_dict


def get_chromosomes_list(seq_dict):
	'''
	From a dictionary, retrieves a list of the keys (the chromosomes).
	'''

	seq_dict_keys = list(seq_dict.keys())
	seq_dict_keys.sort()

	return seq_dict_keys


def process_chromosome(all_reads, offset):
	'''
	Processes through all the reads given, returns a set with Asite sequences.
	'''

	sequence = {}
	for read in all_reads:

		# Reads less than a minimum read length (25) are excluded.

		if read.qlen < 25:
			continue

		protect_nts = read.positions
		protect_nts.sort()

		# support for collapsed reads
		print(int(read.qname.split("_x")[1]))

		if "_x" in read.qname:
			read_count = int(read.qname.split("_x")[1])
		else:
			read_count = 1

		# Aligns to forward strand

		if not read.is_reverse:
			Asite = protect_nts[offset]
		else:
			Asite = protect_nts[-1 - offset] #  -1 for browser display.
		if Asite in sequence:
			sequence[Asite] += read_count
		else:
			sequence[Asite] = read_count
		
	return sequence


def write_bedfile(filepath, seq_dict_keys, offset):
	'''
	Writes a bed file
	'''

	bedfile = open("{0}.bed".format(filepath), "w")
	alignments = pysam.Samfile(filepath, "rb")

	for chrom in seq_dict_keys:
		try:
			all_reads = alignments.fetch(chrom)
		except:
			raise Exception ("Error in fetching ", chrom, " from the alignments")

		sequence = process_chromosome(all_reads, offset)

		for Asite in sorted(sequence):
			bedfile.write("%s\t%s\t%s\t%s\n"%(chrom,  Asite,  Asite+1, sequence[Asite]))
		del sequence
		del all_reads

	bedfile.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=" ")
	parser.add_argument("filepath", type = str, help = "Path to the bam file and relative bai index")
	parser.add_argument("offset", type = str, help = "The offset value. Depends on the treatment and type of experiment")
	parser.add_argument("genomic_fasta", type = str, help = "The genomic fasta")

	args = parser.parse_args()

	seq_dict_keys = get_chromosomes_list(get_fasta_as_dictionary(args.genomic_fasta))

	write_bedfile(args.filepath, seq_dict_keys, int(args.offset))
