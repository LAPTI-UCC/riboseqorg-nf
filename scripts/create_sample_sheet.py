'''
This script takes the RiboSeqOrg Metadata file and a study accession number and creates a sample sheet for the study.

Usage: python create_sample_sheet.py -m <metadata_file> -a <study_accession>

'''

import argparse
import pandas as pd


def main(args):
    metadata = pd.read_csv(args.metadata)
    metadata = metadata[metadata['study_accession'] == args.study_accession]
    metadata = metadata[['study_accession', 'Run', 'ScientificName', 'LIBRARYTYPE']]
    metadata.to_csv(f"{args.output}/sample_sheet.txt", sep='\t', index=False, header=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create sample sheet for a study')
    parser.add_argument('-m', '--metadata', help='RiboSeqOrg Metadata file', required=True)
    parser.add_argument('-s', '--study_accession', help='Study accession number', required=True)
    parser.add_argument('-o', '--output', help='Output file name', required=False)
    args = parser.parse_args()
    main(args)