#!/usr/bin/env python3

import sys
import pysam 
from Bio import SeqIO
import os

filepath = sys.argv[1]
offset = int(sys.argv[2])
genomic_fasta = sys.argv[3]

alignments = pysam.Samfile(filepath,  "rb")
try:
	in_seq_handle = open(genomic_fasta)
except:
	print ("genome file will not open")
	sys.exit()

seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle,  "fasta"))
in_seq_handle.close()
seq_dict_keys = list(seq_dict.keys())
seq_dict_keys.sort()
bedfile = open("{0}.bed".format(filepath), "w")
for chrom in seq_dict_keys:
	chrom_len = len(seq_dict[chrom])
	try:
		all_reads = alignments.fetch(chrom)
	except:
		print (chrom,  "error")
		for i in alignments:
			print(i)
		sys.exit()
	sequence = {}
	for read in all_reads:
		# Reads less than minreadlen are excluded.
		if read.qlen < 25:
			continue
		protect_nts = read.positions
		protect_nts.sort()
		#aligns to forward strand
		if not read.is_reverse:
			Asite = protect_nts[offset]
		else:
			Asite = protect_nts[-1 - offset] #  -1 for browser display.
		if Asite in sequence:
			sequence[Asite] += 1
		else:
			sequence[Asite] = 1
	# print(len(sequence))

	if sequence == {}:
		raise Exception (" There are no valid A sites ")

	for Asite in sorted(sequence):
		bedfile.write("%s\t%s\t%s\t%s\n"%(chrom,  Asite,  Asite+1, sequence[Asite]))
	del sequence
	del all_reads
bedfile.close()
#bedfile = ("{0}{1}.bam_sorted.bam.bed".format(filepath, filename))
#sorted_bedfile = ("{0}{1}.bam_sorted.bam.bed.sorted".format(filepath, filename))
#subprocess.call("sort -k1,1 -k2,2n {0} > {1}".format(bedfile, sorted_bedfile), shell=True)
#subprocess.call("bedtools genomecov -ibam {0}{1}.bam_sorted.bam -g {2} -bg > {0}{1}.cov".format(filepath, filename, master_dict["chrom_sizes_file"]), shell=True)
#subprocess.call("sort -k1,1 -k2,2n {0}{1}.cov > {0}{1}.sorted.cov".format(filepath, filename), shell=True)
