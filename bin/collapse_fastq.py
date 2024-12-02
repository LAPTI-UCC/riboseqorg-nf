import argparse
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def collapse(infile: str, outfile: str) -> None:
    '''
    Read in a fasta file and only store unique sequences
    Output a fasta file with unique sequences and their counts in the header

    Header format: >seq{read_number}_x{count}

    Input:
        infile: str
        outfile: str

    '''
    unique_reads = {}

    with open(infile, 'rb') as f:
        magic_number = f.read(2)
        if magic_number == b'\x1f\x8b':
            f = gzip.open(infile, 'rt')
        else:
            f.seek(0)
            f = open(infile, 'r')

        for title, sequence, quality in FastqGeneralIterator(f):
            if sequence in unique_reads:
                unique_reads[sequence]['count'] += 1
            else:
                unique_reads[sequence] = {}
                unique_reads[sequence]['count'] = 1

    if outfile.endswith('.gz'):
        f = gzip.open(outfile, 'wt')
    else:
        f = open(outfile, 'w')

    read_number = 1
    for seq, vals in unique_reads.items():
        f.write(f'>seq{read_number}_x{vals["count"]}\n')
        f.write(f"{seq}\n")
        read_number += 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='path to input FASTQ file')
    parser.add_argument('-o', help='path to output collapsed FASTQ file (add .gz to the end to compress)')

    args = parser.parse_args()
    collapse(args.i, args.o)