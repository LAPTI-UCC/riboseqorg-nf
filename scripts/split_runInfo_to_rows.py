from cProfile import run
import pandas as pd
import subprocess
import os
import argparse


def run_get_fastq(runInfo_path, outdir):
    '''
    Use ffq to find the url for the fastq and then download to outdir with wget 
    '''
    runInfo = pd.read_csv(runInfo_path, header=0)

    if outdir[-1] != "/":
        outdir = outdir + "/"

    grouped = runInfo.groupby(runInfo['Run'])

    for run in runInfo['Run'].unique():
        temporary_df = grouped.get_group(run)
        temporary_df.to_csv(f"{outdir}/{run}_sraRunInfo.csv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-r", help="Path to sraRunInfo")
    parser.add_argument("-o", help="Path to output directory")


    args = parser.parse_args()
    run_get_fastq(args.r, args.o)
