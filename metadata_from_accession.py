import subprocess
import pandas as pd


def get_geo_results_as_df(accession):
    '''
    Given an accession (GEO SRA or ENA) return the search results from GEO as a pandas data frame. 
    '''
    s = subprocess.check_output(f'esearch -db sra -query {accession} | efetch -format runinfo -mode text', stderr=subprocess.STDOUT, shell=True)
    arr = []
    for i in str(s).split('\\n')[:-1]:
        arr.append(i.split(','))

    arr[0][0] = arr[0][0].strip("b'")

    df = pd.DataFrame(arr[1:], columns=arr[0])
    return(df)


headings = ['Accession', # GEO accession 
    'Title', # Title from GEO submission 
    'Organism', # ; separated list of organisms in this study  
    'Samples', # number of samples in this study 
    'SRA', #SRA project accession (SRP)
    'Release_Date',	# Date of publucation from GEO
    'All protocols', # This column comes from argeos (23-2-22 not yet implemented manually)
    'Type', # Type of sequenicing data according to GEO
    'GSE', # https link to geo entry 
    'GSE_Supplementary', # FTP link to geo entry supplementary data 
    'BioProject', # http link to entry on BioProjetcs 
    'PMID', # Publcations PMID 
    'authors', # Papers authors 
    'abstract', # papers abstract 
    'title', # papers title 
    'doi', # papers doi 
    'date_published', # date the paper was published 
    'PMC', # pmc id 
    'journal' # Journal in which the paper was published
    ]

