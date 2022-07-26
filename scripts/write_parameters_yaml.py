import argparse
import sqlite3
import os
import pandas as pd



def merge_adapter_reports(fastq_dir):
    '''
    Combine all found adapters into one report of the same format.
    '''
    print("Merging adapter reports")

    all_adapters = []

    dir_list = os.listdir(fastq_dir)
    for i in dir_list:
        if i.endswith('adapter_report.tsv'):
            absolute_path = fastq_dir + '/' + i
            with open(absolute_path, 'r') as report:
                lines = report.readlines()
                for line in lines:
                    adapter = line.split('\t')[1].strip('\n')
                    if adapter not in all_adapters: all_adapters.append(adapter)
    
    final_report_path = fastq_dir + '/final_adapter_report.fa'
    with open(final_report_path, 'w') as final_report:
        for idx, adapter in enumerate(all_adapters):
            final_report.write(f">adapter{idx+1}\n{adapter}\n")

    return final_report_path



def get_annotation_organism(path, db='annotation_inventory.sqlite'):
    '''
    Look at the run into file at the given path and find the scientific 
    name of the organism  to which these runs apply.
    Check the annotation inventory for the 'organism' name that applies
    '''
    runInfo = pd.read_csv(path)
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



def write_paramters_yaml(organism, adapter_report_path, yaml_outpath, skip_gwips=False, skip_trips=False, annotations_inventory_sqlite='annotation_inventory.sqlite'):
    '''
    Write paramaters yaml file for given study
    '''
    print('Writing YAML')
    parameter_dict = {'skip_trips':skip_trips, 'skip_gwips':skip_gwips}
    parameter_order = ["adapter_fasta",
                    "rRNA_index",
                    "transcriptome_index",
                    "genome_index",
                    "genome_fasta",
                    "annotation_sqlite",
                    "chrom_sizes_file",
                    "study_dir",
                    "fastq_files",
                    "skip_trips",
                    "skip_gwips"]



    connection = sqlite3.connect("{}".format(annotations_inventory_sqlite))
    cursor = connection.cursor()
    reference_details = cursor.execute(
        f"SELECT * FROM annotation_inventory WHERE organism == '{organism}'"
    ).fetchall()[0]

    parameter_dict['adapter_fasta'] = adapter_report_path
    parameter_dict['study_dir'] = '/'.join(adapter_report_path.split('/')[:-2])
    parameter_dict['fastq_files'] = '/'.join(adapter_report_path.split('/')[:-1]) + '/*.fastq.gz'

    for i in zip(parameter_order[1:-2], reference_details[1:]):
        parameter_dict[i[0]] = i[1]

    with open(yaml_outpath, 'w') as yaml:
        for line in parameter_order:
            if type(parameter_dict[line]) is str:
                yaml.write(f'{line} : "{parameter_dict[line]}"\n')

            elif type(parameter_dict[line]) is bool:
                yaml.write(f'{line} : {str(parameter_dict[line]).lower()}\n')


def run_project_setup(adapter_report_dir, sraRunInfo, paramters_yaml_path, db='annotation_inventory.sqlite'):
    '''
    given a sra run info file for a study, run the necessary data fetching and setup. 
    '''
    final_adapter_report_path = merge_adapter_reports(adapter_report_dir)
    annotation_inventory_organism = get_annotation_organism(sraRunInfo, db)
    write_paramters_yaml(annotation_inventory_organism, final_adapter_report_path, paramters_yaml_path, annotations_inventory_sqlite=db)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", help="Path to the adapter report directory")
    parser.add_argument("-s", help="path to sqlite of annotation inventory", default="annotation_inventory.sqlite")
    parser.add_argument("-r", help="runInfo file for this study from SRA")
    parser.add_argument("-o", help="outpath for paramters.yaml")

    args = parser.parse_args()
    run_project_setup(args.a, args.r, args.o, db=args.s)