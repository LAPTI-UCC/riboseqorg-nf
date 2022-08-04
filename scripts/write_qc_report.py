
import argparse
import sqlite3
from sqlitedict import SqliteDict
import pickle5
import sys 
import os


def parse_fastqc(report_path):
    '''
    Open the fastqc and return a dictionary of report records {key:string}
    '''
    parsed_fastqc = {}
    with open(report_path, 'r') as file:
        lines = file.readlines()
        
        for line in lines:
            if line.startswith('>>END_MODULE'):
                parsed_fastqc[module] = module_content

            elif line.startswith('>>'):
                module_line = line.strip('>>').split('\t')
                module = module_line[0]
                module_content = module_line[1]

            elif line.split('\t')[0] != '##FastQC':
                module_content += line

    return parsed_fastqc



def process_fastqc(report_path):
    '''
    Obtain desired metrics from fastqc reprot at provided path
    '''
    qc_report = {}

    parsed_fastqc = parse_fastqc(report_path)

    # for entry in parsed_fastqc:
    #     print(entry)
    #     for line in parsed_fastqc[entry].split('\n'):
    #         print(line)

    return parsed_fastqc


def my_decoder(obj):
	return pickle5.loads(obj)


def read_sqlite_dict(filepath):
    '''
    read in sqlite dict from filepath
    '''
    sqlite_db = SqliteDict(f"{filepath}", autocommit=False, decode=my_decoder)
    return sqlite_db


def calculate_triplet_periodicity_score(sqlite_dict, min_read_length=25, max_read_length=35):
    '''
    Triplet periodicity score for riboseq study
    1-(count of the second highest peak/count of the highest peak), 
    if more than one readlength is shown the highest peaks and second highest 
    peaks of each readlength is summed and then the formula is applied.
    '''
    highest, second_highest = 0, 0
    for readlength in range(min_read_length, max_read_length+1):
        sorted_frame_counts = sorted([sqlite_dict["fiveprime"][readlength]['0'], 
                                        sqlite_dict["fiveprime"][readlength]['1'], 
                                        sqlite_dict["fiveprime"][readlength]['2']])
        highest += sorted_frame_counts[2]
        second_highest += sorted_frame_counts[1]

    trip_periodicity_score = 1 - (second_highest/highest)
    return round(trip_periodicity_score , 4)


def get_sqlite_cursor(organism_sqlite_path):
    '''
    get a cursor for the sqlite at the provided path 
    '''
    conn = sqlite3.connect(organism_sqlite_path)
    return conn.cursor()


def query_database(query, cursor):
    '''
    return the results of the provided query
    '''
    results = cursor.execute(query).fetchall()
    return results



def get_average_readlengths(read_length_dict):
    '''
    from read_lengths dict return the average read length 
    '''
    numerator = 0
    number_of_reads = 0
    for key in read_length_dict:
        number_of_reads += read_length_dict[key]
        numerator += key * read_length_dict[key]
    return round(numerator/number_of_reads, 2)


def process_transcript_reads(transcript, cds_start, cds_stop, sequence_length, transcript_reads, offsets):
    '''
    get gene body distribution and profile for this transcript
    '''
    gene_body = {'fiveprime':0, 'cds':0, 'threeprime':0}

    metagene = {'fiveprime':{i:0 for i in range(-300, 301)}, 
                    'threeprime':{i:0 for i in range(-300, 301)}}

    for readlength in transcript_reads:
        for position in transcript_reads[readlength]:
            adjusted_position = position + offsets['fiveprime']["offsets"][readlength]
            if adjusted_position in range(0, cds_start):
                
                gene_body['fiveprime'] += transcript_reads[readlength][position]

            elif adjusted_position in range(cds_start, cds_stop):
                gene_body['cds'] += transcript_reads[readlength][position]

            elif adjusted_position in range(cds_stop, sequence_length):
                gene_body['threeprime'] += transcript_reads[readlength][position]
            
            if adjusted_position in range(cds_start-300, cds_start+300):
                relative_position = adjusted_position - cds_start
                metagene['fiveprime'][relative_position] += transcript_reads[readlength][position]
            
            elif adjusted_position in range(cds_stop-300, cds_stop+300):
                relative_position = adjusted_position - cds_stop
                metagene['threeprime'][relative_position] += transcript_reads[readlength][position]

    return gene_body, metagene


def generate_profile(sqlite_dict, organism_sqlite):
    '''
    create a whole gene metagene profile from the given sqlite dict
    reads too far from start or stop are combined in one variable.
    '''
    cursor = get_sqlite_cursor(organism_sqlite)
    transcript_table = query_database("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE principal = 1;", cursor)
    offsets = sqlite_dict['offsets']
    gene_body = {'fiveprime':0, 'cds':0, 'threeprime':0}

    metagene = {'fiveprime':{i:0 for i in range(-300, 301)}, 
                'threeprime':{i:0 for i in range(-300, 301)}}

    counter = 0
    for row in transcript_table:
        transcript, cds_start, cds_stop, sequence = row
        if transcript in sqlite_dict:
            transcript_reads = sqlite_dict[transcript]["unambig"]
        else:
            continue

        if cds_start != None and cds_stop !=None:
            transcript_gene_body, transcript_metagene = process_transcript_reads(transcript, cds_start, cds_stop, len(sequence), transcript_reads, offsets)
        
            for key in transcript_gene_body:
                gene_body[key] += transcript_gene_body[key]

            for key in transcript_metagene:
                for position in transcript_metagene[key]:
                    metagene[key][position] += transcript_metagene[key][position]


        counter += 1
    return metagene, gene_body


def sliding_window(elements, window_size=100):
    '''
    given a list of elements and a set window size return a nested list of sliding windows 
    '''
    if len(elements) <= window_size:
       return None

    window_list = []
    for i in range(len(elements)):
        window_list.append(elements[i:i+window_size])
    
    return window_list


def calculate_metagene_summary_score(metagene):
    '''
    Generate a summary score for both metagene profiles in metagene dict
    '''

    for prime in metagene:
        print()
        print('fiveprime', sum([i for i in metagene[prime] if i < 0])/300)
        print('threeprime', sum([i for i in metagene[prime] if i >= 0])/300)


        positions = sorted(metagene[prime].keys())
        for window in sliding_window(positions):
            metagene_window = [metagene[prime][i] for i in window]
            # print(sum(metagene_window)/len(metagene_window))



def process_readfile(readfile_path, organism_sqlite):
    '''
    Get desired metrics from read file 
    '''
    readfile_report = {}

    sqlite_dict = read_sqlite_dict(readfile_path)

    if "trip_periodicity" in sqlite_dict:
        trip_peridicity = sqlite_dict['trip_periodicity']
        trip_periodicity_score = calculate_triplet_periodicity_score(trip_peridicity)
    else:
        raise Exception(f"No triplet periodicity in sqlite database {readfile_path}")

    if "read_lengths" in sqlite_dict:
        read_lengths = sqlite_dict["read_lengths"]
        average_readlength = get_average_readlengths(read_lengths)
    else:
        raise Exception(f"No read length distribution data in sqlite database {readfile_path}")

    metagene, gene_body = generate_profile(sqlite_dict, organism_sqlite)
    
    calculate_metagene_summary_score(metagene)

    gene_body_proportion_cds = round(gene_body['cds']/(gene_body['fiveprime'] + gene_body['cds'] + gene_body['threeprime']), 2)
    gene_body_ratio = round(gene_body['cds']/(gene_body['fiveprime'] + gene_body['threeprime']), 2)

    return readfile_report


def write_final(qc_report, readfile_report, outpath):
    '''
    combine the metrics to a standard report. 
    '''

    with open(outpath, 'w') as outfile:
        for entry in qc_report:
            outfile.write(f"{entry} : {qc_report[entry]}")

        for entry in readfile_report:
            outfile.write(f"{entry} : {readfile_report[entry]}")



def main(report_path, readfile_path, organism_sqlite, outpath):
    '''
    produce qc report for this dataset and write it to outpath
    '''
    qc_report = process_fastqc(report_path)

    readfile_report = process_readfile(readfile_path, organism_sqlite)

    write_final(qc_report, readfile_report, outpath)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", help="path to sqlite readfile")
    parser.add_argument("-r", help="path to fastqc report")
    parser.add_argument("-p", help="path to sqlite for organism")
    parser.add_argument("-o", help="output path for the report")

    args = parser.parse_args()
    main(args.r, args.s, args.p, args.o)


'''

The output of this script shoud be like an extended fastqc report with added ribo-seq metrics

Include fastqc report at top of file and then a similar style riboqc report

- basic statistics (averages or scores for each stat)
- triplet periodicity (done)
- read length distribution (done)
- gene body distribution 
- metagene start counts
- metagene stop counts 

'''