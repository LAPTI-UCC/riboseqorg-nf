'''
This script takes the RiboSeqOrg Metadata file and a study accession number and creates a sample sheet for the study.

Usage: python create_sample_sheet.py -m <metadata_file> -a <study_accession>

'''

import argparse
import pandas as pd

import os
import glob


def get_processed_projects(path="/data1/Jack/projects/Collapse-FASTQ/data/RiboSeqOrg/ribometric/"):
    '''
    Get a list of projects that have been processed
    
    Keyword Arguments:
        path {str} -- Path to the processed data (default: {"/data1/Jack/projects/Collapse-FASTQ/data/RiboSeqOrg/ribometric/"})
    
    Returns:
        [list] -- List of projects that have been processed
    '''
    projects = os.listdir(path)
    projects = [x.split('.')[0] for x in projects]
    return projects


def get_fasta_path(
        run: str,
        bioproject: str,
        base_path="/data1/Jack/projects/Collapse-FASTQ/data/RiboSeqOrg/collapsed_fastq/",
        ext="_clipped_collapsed.fastq.gz",):
    '''
    Using OS find a valid path for each run
    
    Arguments:
        run {str} -- Run accession
        bioproject {str} -- BioProject accession

    Keyword Arguments:
        base_path {str} -- Base path to the data (default: {"/data1/Jack/projects/Collapse-FASTQ/data/RiboSeqOrg/collapsed_fastq/"})
        ext {str} -- Extension of the file (default: {"_clipped_collapsed.fastq.gz "})

    Returns:
        [str] -- Valid path to the file
    '''
    matching_files = glob.glob(os.path.join(base_path, '**', f'{run}*{ext}'), recursive=True)
    return matching_files[0] if len(matching_files) > 0 else None


def main(args):
    metadata = pd.read_csv(args.metadata)

    if args.species:
        metadata = metadata[metadata['ScientificName'] == args.species]
    
    metadata = metadata[metadata['LIBRARYTYPE'] == 'RFP']

    processed_projects = get_processed_projects()
    metadata = metadata[~metadata['BioProject'].isin(processed_projects)]


    # # get project ids  
    bioprojects = metadata['BioProject'].unique()[:args.num_projects]
    metadata = metadata[metadata['BioProject'].isin(bioprojects)]

    metadata['path'] = metadata.apply(lambda x: get_fasta_path(x['Run'], x['BioProject']), axis=1)

    metadata = metadata[[ 'Run', 'BioProject', 'ScientificName', 'path']]

    # print number of Nones in path

    metadata.to_csv(f"/data1/riboseq_org/riboseq_data_processing/data/RiboSeqOrg/sample_sheet.csv", sep='\t', index=False, header=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create sample sheet for a study')
    parser.add_argument('-m', '--metadata', help='RiboSeqOrg Metadata file', required=True)
    parser.add_argument('-s', '--species', help='Species name', required=False)
    parser.add_argument('-o', '--output', help='Output dir name', required=False)
    parser.add_argument('-n', '--num_projects' , help='Number of projects to process', required=False, default=10, type=int)
    args = parser.parse_args()
    main(args)