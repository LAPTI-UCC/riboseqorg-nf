import argparse
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def collapse(infile, outfile):
    # Store the unique reads in a dictionary with their counts
    unique_reads = {}

    # Open the FASTQ file for reading, handling gzip if necessary
    # by checking the magic number at the beginning of the file
    with open(infile, 'rb') as f:
        magic_number = f.read(2)
        if magic_number == b'\x1f\x8b':
            f = gzip.open(infile, 'rt')
        else:
            f.seek(0)
            f = open(infile, 'r')

        # Read the file line by line
        for title, sequence, quality in FastqGeneralIterator(f):
            if sequence in unique_reads:
                unique_reads[sequence]['count'] += 1
            else:
                unique_reads[sequence] = {}
                unique_reads[sequence]['count'] = 1
                unique_reads[sequence]['qual'] = quality


        # for line in map(str.strip, f):
        #     # If the line is a sequence, add it to the dictionary of unique reads
        #     # and increment its count
        #     if line.startswith('@'):
        #         seq = next(f)
        #         if seq in unique_reads:
        #             # unique_reads[seq] += 1
        #             unique_reads[seq]['count'] += 1
        #             # unique_reads[seq]['qual']
        #             next(f)

        #         else:
        #             next(f)
        #             # unique_reads[seq] = 1
        #             unique_reads[seq] = {}
        #             unique_reads[seq]['count'] = 1
        #             unique_reads[seq]['qual'] = next(f)
        #         # Skip the remaining lines for this read
        #         next(f)
        #         # print(unique_reads[seq])
    # Open the FASTQ file for writing, handling gzip if necessary
    # by checking the magic number at the beginning of the file
    if outfile.endswith('.gz'):
        f = gzip.open(outfile, 'wt')
    else:
        f = open(outfile, 'w')

        # Write the unique reads to the file with their counts
    read_number = 1
    for seq, vals in unique_reads.items():
        f.write(f'@read{read_number}_x{vals["count"]}\n')
        f.write(f"{seq}\n")
        f.write('+\n')
        f.write(f'{vals["qual"]}\n')
        read_number += 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='path to input FASTQ file')
    parser.add_argument('-o', help='path to output collapsed FASTQ file (add .gz to the end to compress)')

    args = parser.parse_args()
    collapse(args.i, args.o)