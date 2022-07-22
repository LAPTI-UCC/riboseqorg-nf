import subprocess
import argparse
import pandas as pd
import os
import json

from simplejson import loads 


def run_get_fastq(runInfo_path, outdir):
    '''
    Use ffq to find the url for the fastq and then download to outdir with wget 
    '''
    runInfo = pd.read_csv(runInfo_path, header=0)

    if outdir[-1] != "/":
        outdir = outdir + "/"

    for idx, row in runInfo.iterrows():
        ffq_stdout = subprocess.run(f"ffq --ftp {row['Run']} | jq -r .[]", check=True, capture_output=True, shell=True)

        ffq_metadata_dict = json.loads(ffq_stdout.stdout.decode())

        # if not os.path.exists(outdir):
        #     os.makedirs(outdir)
        a = subprocess.run(f"wget {ffq_metadata_dict['url']} -P {outdir} ", check=True, capture_output=True, shell=True)
        return a

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", help="runInfo file for this study from SRA")
    parser.add_argument("-o", help="outpath (Dir) for resulting fastqs")


    args = parser.parse_args()
    run_get_fastq(args.r, args.o)