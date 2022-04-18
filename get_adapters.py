import sys
import subprocess



def rev_comp(adapter):
    '''
    Return the reverse complimnet of a provided nucleotide adapter sequence 
    '''
    adapter = adapter[::-1]
    adapter = (
        adapter.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g")
    )
    return adapter.upper()


def check_adapter(adapter, fastq_path, verbose=False):
    '''
    For a given adapter sequence check for its presence in the given FASTQ file 
    '''
    fastq_lines_output = subprocess.check_output(
        "wc -l {0}".format(fastq_path), shell=True
    )
    fastq_lines = float(fastq_lines_output.split()[0])
    number_of_reads = fastq_lines/4


    adapter_count_raw = subprocess.check_output(
        "head -2000000 {0} | sed -n '2~4p' > test.fastq ; agrep -c1 \"{1}\" test.fastq ; rm test.fastq".format(
            fastq_path, adapter
        ),
        shell=True,
    )

    adapter_count = float(adapter_count_raw.decode('utf-8').strip('\n'))

    percentage_contamination = float((adapter_count / number_of_reads) * 100)
    if percentage_contamination >= (0.05):
        return True
    else:
        return False


def get_adapters(fastq_path, adapter_sequences, verbose=False):
    '''
    Given a list of know adapters check if each (or its reverse compliment) is found in the fastq file for which the path was provided
    '''

    found_adapters ={'forward': [], 'reverse': []}

    for adapter in adapter_sequences:
        verdict = check_adapter(adapter, fastq_path, verbose=verbose)
        if verdict:
            found_adapters['forward'].append(adapter)

        else:
            adapter = rev_comp(adapter)
            verdict = check_adapter(adapter, fastq_path, verbose=verbose)
            if verdict:
                found_adapters['reverse'].append(adapter)

    return found_adapters

    

def write_adapter_report(found_adapters, outfile_path):
    '''
    write a simple output file that lists the found forward and reverse adapters 
    '''
    with open(outfile_path, 'w') as outfile:

        for direction in found_adapters:
            for adapter in found_adapters[direction]:
                outfile.write(f'{direction}\t{adapter}\n')




if __name__ == '__main__':
    fastq_path = sys.argv[1]
    path_list = fastq_path.split('/')
    fastq_dir = '/'.join(path_list[:-1])
    fastq_filename = path_list[-1]

    adapter_sequences = [
    "CTGTAGGCACCATCAAT",
    "AGATCGGAAGAGC",
    "CGCCTTGGCCGTACAGCAG",
    "AAAAAAAAAAAAA",
    "TGGAATTCTCGGGTGCCAAGG",
    "CCTTGGCACCCGAGAATT",
    "GATCGGAAGAGCGTCGT",
    "CTGATGGCGCGAGGGAG",
    "GATCGGAAGAGCACACG",
    "AATGATACGGCGACCAC",
    "GATCGGAAGAGCTCGTA",
    "CAAGCAGAAGACGGCAT",
    "ACACTCTTTCCCTACA",
    "GATCGGAAGAGCGGTT",
    "ACAGGTTCAGAGTTCTA",
    "CAAGCAGAAGACGGCAT",
    "ACAGGTTCAGAGTTCTA",
    "CAAGCAGAAGACGGCAT",
    "TGATCGGAAGAGCACAC",
    ]


    found_adapters = get_adapters(fastq_path=fastq_path, adapter_sequences=adapter_sequences)

    report_path = f'{fastq_dir}/{fastq_filename}_adapter_report.tsv'
    write_adapter_report(found_adapters, report_path)
