from logging import raiseExceptions
import pandas as pd
import subprocess
import os
import sqlite3
import argparse
from copy import copy



def get_study_structure(runInfo_path):
    '''
    From the organism names in the run into column ScientificNames create subdirectories
    for each sub-study (ie. a study that uses just one organism.)
    Also return the list of paths to the RunInfos that are put in these directories 
    '''
    runInfo_paths = []

    study_dir = os.path.dirname(runInfo_path)

    runInfo_file_name = os.path.basename(runInfo_path)
    runInfo_file_name_parts = runInfo_file_name.split('_')
    keep = [runInfo_file_name_parts[0]] 
    for part in runInfo_file_name_parts:
        if 'GSE' in part or 'SRP' in part: 
            keep.append(part)

    keep.append(runInfo_file_name_parts[-1])

    runInfo = pd.read_csv(runInfo_path)
    if runInfo['ScientificName'].nunique() == 1:
        runInfo_paths.append(runInfo_path)

    else:
        for i in runInfo['ScientificName'].unique():
            new_runInfo = runInfo[runInfo['ScientificName'] == i]
            organism_name = '_'.join(i.split(' '))

            new_runInfo_namelist = copy(keep)
            new_runInfo_namelist.insert(1, organism_name)
            new_runInfo_filename = '_'.join(new_runInfo_namelist)
            new_study_dir = study_dir + '/' + new_runInfo_filename.strip('_sraRunInfo.csv')
            if not os.path.exists(new_study_dir):
                os.makedirs(new_study_dir)

            new_runInfo.to_csv(f"{new_study_dir}/{new_runInfo_filename}" )
            runInfo_paths.append(f"{new_study_dir}/{new_runInfo_filename}")

    return runInfo_paths

def download_files_from_SRA(runInfo_path, outdir):
    '''
    download all runs to a specifc output dir
    '''
    print("Downloading from SRA")
    runInfo = pd.read_csv(runInfo_path, header=0)

    if outdir[-1] != "/":
        outdir = outdir + "/"

    for idx, row in runInfo.iterrows():
        if not os.path.isfile(outdir + row['Run']):
            subprocess.run(['wget', row['download_path'], '-P', outdir])


def sra_to_fastq(study_dir):
    '''
    convert all .sra files in sra_dir to fastq
    '''
    print("Converting SRA to FQ")
    sra_dir = study_dir + '/sra'
    dir_list = os.listdir(sra_dir)
    for i in dir_list:
        absolute_path = sra_dir + '/' + i
        
        if not os.path.isfile(study_dir + '/fastq/' + i + '.fastq'):
            subprocess.run(['fasterq-dump', absolute_path, '-O', study_dir+'/fastq'])


def find_adapters(fastq_dir):
    '''
    check fastq files to see what adapters are present.
    '''
    print("Finding Adapters")
    dir_list = os.listdir(fastq_dir)
    for i in dir_list:
        if i.endswith('.fastq'):
            absolute_path = fastq_dir + '/' + i
            if not os.path.isfile(absolute_path + '_adapter_report.tsv'):
                print("Finding Adapters: ", absolute_path)
                subprocess.run(['python', '/home/jack/projects/riboseq_data_processing/get_adapters.py', absolute_path])


def merge_adapter_reports(fastq_dir):
    '''
    Combine all found adapters into one report of the same format.
    '''
    print("Merging adapter reports")

    all_adapters = []

    dir_list = os.listdir(fastq_dir)
    for i in dir_list:
        if i.endswith('.fastq_adapter_report.tsv'):
            absolute_path = fastq_dir + '/' + i
            with open(absolute_path, 'r') as report:
                lines = report.readlines()
                for line in lines:
                    adapter = line.split('\t')[1].strip('\n')
                    if adapter not in all_adapters: all_adapters.append(adapter)
    
    with open(fastq_dir + '/final_adapter_report.tsv', 'w') as final_report:
        for idx, adapter in enumerate(all_adapters):
            final_report.write(f"adapter{idx+1}\t{adapter}\n")


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
    print('Writing XML')
    parameter_dict = {'skip_trips':skip_trips, 'skip_gwips':skip_gwips}
    parameter_order = ["adapter1",
                    "adapter2",
                    "rRNA_index",
                    "transcriptome_index",
                    "genome_index",
                    "genome_fasta",
                    "annotation_sqlite",
                    "chrom_sizes_file",
                    "skip_trips",
                    "skip_gwips"]

    with open(adapter_report_path, 'r') as adapter_report:
        lines = adapter_report.readlines()
        if len(lines) == 0:
            parameter_dict["adapter1"], parameter_dict["adapter2"] = "", ""
        elif len(lines) <= 2:
            parameter_dict["adapter1"] = lines[0].split("\t")[1].strip('\n')
            parameter_dict["adapter2"] = lines[1].split("\t")[1].strip('\n')


    connection = sqlite3.connect("{}".format(annotations_inventory_sqlite))
    cursor = connection.cursor()
    reference_details = cursor.execute(
        f"SELECT * FROM annotation_inventory WHERE organism == '{organism}'"
    ).fetchall()[0]

    for i in zip(parameter_order[2:-2], reference_details[1:]):
        parameter_dict[i[0]] = i[1]

    with open(yaml_outpath, 'w') as yaml:
        for line in parameter_order:
            if type(parameter_dict[line]) is str:
                yaml.write(f'{line} \t:\t "{parameter_dict[line]}"\n')

            elif type(parameter_dict[line]) is bool:
                yaml.write(f'{line} \t:\t {str(parameter_dict[line]).lower()}\n')



def run_project_setup(run_info, db='annotation_inventory.sqlite'):
    '''
    given a sra run info file for a study, run the necessary data fetching and setup. 
    '''

    runInfo_paths = get_study_structure(run_info)
    for path in runInfo_paths:
        print('-'*90, '\n')

        study_dir = os.path.dirname(path)

        download_files_from_SRA(path, f'{study_dir}/sra')
        sra_to_fastq(study_dir)

        find_adapters(f'{study_dir}/fastq')
        merge_adapter_reports(study_dir+'/fastq')
        try:
            annotation_inventory_organism = get_annotation_organism(path, db='annotation_inventory.sqlite')
            write_paramters_yaml(annotation_inventory_organism, study_dir + '/fastq/final_adapter_report.tsv', study_dir + '/parameters.yaml')
        except:
            print(annotation_inventory_organism)
            continue


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", help="The organism from which the samples have come")
    parser.add_argument("-s", help="path to sqlite of annotation inventory", default="annotation_inventory.sqlite")
    parser.add_argument("-r", help="runInfo file for this study from SRA")

    args = parser.parse_args()
    run_project_setup(args.r, db=args.s)

