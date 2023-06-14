'''
Given a sample sheet, create a params file for the pipeline

Usage: python create_params_file.py -s <sample_sheet> -o <output_file> -a <annotation_inventory>

'''

import argparse
import pandas as pd
import sqlite3


def get_annotation_organism(path, db='annotation_inventory.sqlite'):
    '''
    Look at the run into file at the given path and find the scientific 
    name of the organism  to which these runs apply.
    Check the annotation inventory for the 'organism' name that applies
    '''
    runInfo = pd.read_csv(path, sep='\t')
    scientific_name = runInfo['ScientificName'].unique()[0]

    connection = sqlite3.connect(db)
    cursor = connection.cursor()
    result = cursor.execute(f"SELECT organism FROM primary_organism WHERE scientific_name = '{scientific_name}';").fetchall()
    
    if len(result) > 0:
        return result[0][0]
    elif len(result) == 0:
        return Exception(f"""\tERROR: There is no primary organism assigned for {scientific_name}.\n
         Add one using annnotation inventory managers --operation set_primary_organism
         """)



def fetch_studywide_info(organism,
                         skip_gwips=False,
                         skip_trips=False, 
                         annotations_inventory_sqlite='annotation_inventory.sqlite'
                         ) -> dict:
    '''
    Fetch study wide information from annotation inventory
    this includes the paths to the reference files etc 

    Input:
        organism: str
            The organism name
        skip_gwips: bool
            Whether to skip gwips in the pipeline
        skip_trips: bool
            Whether to skip trips in the pipeline
        annotations_inventory_sqlite: str
            Path to the annotation inventory sqlite database
    '''
    parameter_dict = {'skip_trips':skip_trips, 'skip_gwips':skip_gwips}
    parameter_order = [
                    "rRNA_index",
                    "transcriptome_index",
                    "genome_index",
                    "genome_fasta",
                    "annotation_sqlite",
                    "chrom_sizes_file",
                    "study_dir",
                    "skip_trips",
                    "skip_gwips"
                    ]

    connection = sqlite3.connect("{}".format(annotations_inventory_sqlite))
    cursor = connection.cursor()
    print(organism)
    reference_details = cursor.execute(
        f"SELECT * FROM annotation_inventory WHERE organism == '{organism}'"
    ).fetchall()[0]
    connection.close()
    parameter_dict.update(dict(zip(parameter_order[:-2], reference_details[1:])))
    return parameter_dict


def main(args):
    organism = get_annotation_organism(args.sample_sheet, args.annotation_inventory)
    parameter_dict = fetch_studywide_info(organism,
                                          skip_gwips=args.skip_gwips,
                                          skip_trips=args.skip_trips,
                                          annotations_inventory_sqlite=args.annotation_inventory)
    parameter_dict['study_dir'] = '/'.join(args.sample_sheet.split('/')[:-1])
    parameter_dict['sample_sheet'] = args.sample_sheet

    with open(f"{args.output_dir}/params.config", 'w') as f:
        f.write("params {\n")
        f.write(f'\tstudy_dir = "{parameter_dict["study_dir"]}"\n')
        f.write(f'\tskip_trips = {str(parameter_dict["skip_trips"]).lower()}\n')
        f.write(f'\tskip_gwips = {str(parameter_dict["skip_gwips"]).lower()}\n')
        f.write(f'\ttranscriptome_index = "{parameter_dict["transcriptome_index"]}"\n')
        f.write(f'\tgenome_index = "{parameter_dict["genome_index"]}"\n')
        f.write(f'\tgenome_fasta = "{parameter_dict["genome_fasta"]}"\n')
        f.write(f'\tannotation_sqlite = "{parameter_dict["annotation_sqlite"]}"\n')
        f.write(f'\tchrom_sizes_file = "{parameter_dict["chrom_sizes_file"]}"\n')
        f.write(f'\trRNA_index = "{parameter_dict["rRNA_index"]}"\n')
        f.write(f'\tsample_sheet = "{parameter_dict["sample_sheet"]}"\n')
        f.write("}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create params file for a study')
    parser.add_argument('-s', '--sample_sheet', help='Sample sheet', required=True)
    parser.add_argument('-a', '--annotation_inventory', help='Path to Annotation Inventory', required=False)
    parser.add_argument('-o', '--output_dir', help='Output dir name', required=True)
    parser.add_argument('-g', '--skip_gwips', help='Skip gwips', action='store_true')
    parser.add_argument('-t', '--skip_trips', help='Skip trips', action='store_true')
    args = parser.parse_args()
    main(args)