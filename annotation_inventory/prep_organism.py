'''
Prepare annotation files for a given organism for addition to the annotation inventory. 

Uses the Ensembl REST API to obtain files for the given organism.
'''

import argparse
import requests
import sys
import subprocess

def fetch_endpoint(server: str, request: str, content_type: str ):
    """
    Fetch an endpoint from the server, allow overriding of default content-type
    """
    r = requests.get(server+request, headers={ "Accept" : content_type})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text
    
    
def fetch_endpoint_POST(server: str, request: str, data, content_type='application/json'):

    r = requests.post(server+request,
                      headers={ "Content-Type" : content_type},
                      data=data )

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text
    
    
def get_organism_info(organism: str):
    '''
    Get the Ensembl ID for the given organism.
    
    Args:
        organism (str): The organism to get the Ensembl ID for.
        
    Returns:
        str: The Ensembl ID for the given organism.
    '''
    
    # get the Ensembl ID for the given organism
    server = "http://rest.ensembl.org"
    ext = "/info/assembly/"
    r = requests.get(server+ext+organism, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    
    return decoded['assembly_accession']

def get_organism_gtf(organism: str, output_dir: str):
    '''
    Get the GTF file for the given organism.
    
    Args:
        organism (str): The organism to get the GTF file for.
        output_dir (str): The directory to write the GTF file to.
        
    Returns:
        str: The path to the GTF file.
    '''
    
    # get the GTF file for the given organism
    server = "http://rest.ensembl.org"
    ext = "/sequence/id/"
    r = requests.get(server+ext+organism, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()    
    return decoded['species']

def main(args):
    '''
    Wrapper function for organism prep
    '''
    accession = get_organism_info(args.o)

    # use subprocess to run the following command
    # datasets download genome accession {accession}
    
    subprocess.run(['datasets', 'download', 'genome', 'accession', accession, '-o', args.d])





    return 0

if __name__ == '__main__':
    # set up argparse with params for organism and output directory
    parser = argparse.ArgumentParser(description='Prepare annotation files for a given organism for addition to the annotation inventory.')
    parser.add_argument('-o', help='The organism to prepare annotation files for.')
    parser.add_argument('-d', help='The directory to write the annotation files to.')
    args = parser.parse_args()

    main(args)