import pandas as pd
import subprocess
import os
import argparse

from csv_header_check import check_and_assign_header_from_file, csv_header_check


def get_sra_run_info_as_df(srp, outfile_path):
    """
    Given a srp accession return the sra run table as a pandas data frame.
    """
    subprocess.check_output(
        f"esearch -db sra -query {srp} | efetch -format runinfo -mode text | cat > {outfile_path}",
        stderr=subprocess.STDOUT,
        shell=True,
    )



def run_table_to_csv(sra_runs, outfile):
    """
    Given a pandas df of SRA run table. Write it to csv
    """
    sra_runs.to_csv(outfile)


def superset_row_to_outfile_name(row):
    """
    take a row of the ribosme profiling superset and return a string with a useful file name
    """

    if str(row["authors"]) != "nan":
        first_author_last_name = row["authors"].split(",")[0].split(" ")[-1]
    else:
        first_author_last_name = ""

    if str(row["Release_Date"]) != "nan":
        year = str(row["Release_Date"]).split("-")[0]
    else:
        year = ""

    if ";" not in row["Organism"]:
        outfile = (
            first_author_last_name
            + year
            + "_"
            + "_".join(row["Organism"].split(" "))
            + "_"
            + row["Accession"]
            + "_"
            + row["SRA"]
        )
    else:
        organism_names = []
        organisms = row["Organism"].split(";")
        for i in range(len(organisms)):
            organism_names.append("_".join(organisms[i].split(" ")))

        outfile = (
            first_author_last_name
            + year
            + "_"
            + "".join(organism_names)
            + "_"
            + row["Accession"]
            + "_"
            + row["SRA"]
        )

    return outfile


def download_sra_run_table(
    superset_path, datadir, outfile_path, num_samples=1, specific_GSEs=[], sep=","
):
    """
    download the sra run table for either a random srp or a specified list of srps
    """
    superset = pd.read_csv(superset_path, sep=sep, header=0)
    if datadir[-1] == "/":
        datadir = datadir.strip('/')

    if not specific_GSEs:
        srps = superset.sample(num_samples)
    else:
        srps = superset[superset.Accession.isin(specific_GSEs)]

    for idx, row in srps.iterrows():
        outfile = superset_row_to_outfile_name(row)
        if not os.path.isdir(datadir + "/" + outfile):
            os.makedirs(datadir + "/" + outfile)
        # outfile_path = f"{datadir}/{outfile}/{outfile}_sraRunInfo.csv"
        get_sra_run_info_as_df(row["SRA"], outfile_path)
        print(f"File outputted to {outfile_path}")

        check_and_assign_header_from_file(outfile_path)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "path_to_superset",
        type=str,
        help="The path to the ribosome_profiling_superset.csv",
    )
    parser.add_argument(
        "path_to_datadir", type=str, help="The folder where the data is"
    )
    parser.add_argument(
        "gse_to_fetch", type=str, help="The specific GSE to fetch (es.GSE97384)"
    )
    parser.add_argument(
        "outpath", type=str, help="relative path to output"
    )

    args = parser.parse_args()
    
    download_sra_run_table(
        args.path_to_superset, args.path_to_datadir, args.outpath, specific_GSEs=[args.gse_to_fetch]
    )

    # download_sra_run_table('data/ribosome_profiling_superset.csv',
    #         'data',
    #         specific_GSEs=['GSE97384']
    #         )
