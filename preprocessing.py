import pandas as pd
import subprocess
import os
import sqlite3
import argparse

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
        if len(lines) <= 2:
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
                yaml.write(f'{line} \t:\t "{parameter_dict[line]}"')

            elif type(parameter_dict[line]) is bool:
                yaml.write(f'{line} \t:\t {str(parameter_dict[line]).lower()}')



def run_project_setup(run_info):
    '''
    given a sra run info file for a study, run the necessary data fetching and setup. 
    '''
    study_dir = os.path.dirname(run_info)
    download_files_from_SRA(run_info, f'{study_dir}/sra')
    sra_to_fastq(study_dir)

    find_adapters(f'{study_dir}/fastq')
    merge_adapter_reports(study_dir+'/fastq')

    '''
    The remaining task here is to find the organisms used in the run info and to find the organism name in the database. ]]
    I will also need to split studies with multiple organisms into different study dirs
    '''
    write_paramters_yaml('human_hg19_gencode25', study_dir + '/fastq/final_adapter_report', study_dir + 'parameters.yaml')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", help="The organism from which the samples have come")
    parser.add_argument("-s", help="path to sqlite of annotation inventory", default="annotation_inventory.sqlite")
    parser.add_argument("-r", help="runInfo file for this study from SRA")

    args = parser.parse_args()
    print(args.r)
    run_project_setup(args.r)

