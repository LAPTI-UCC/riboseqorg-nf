import pandas as pd
import argparse

def get_GSE_not_downloaded(df):
    '''
    Returns the list of GSEs the process GET_GSE_REPORT was not able to download
    '''
    failed_df = df.loc[df['status'] == "FAILED"]
    df2 = failed_df[failed_df['name'].str.contains('REPORT')]

    not_downloaded_list = []
    for n in df2.index:
        raw_GSE = df2.at[n,"raw_GSE"]
        GSE = raw_GSE[raw_GSE.index("GSE"):raw_GSE.index("_family")]
        not_downloaded_list.append(GSE)

    return not_downloaded_list

# Need to write this output to a txt file?

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieves the GSEs that have not been downloaded correctly")
    parser.add_argument("path_to_run_report", type = str, help = "Path to the location of the run report, a csv file")
    args = parser.parse_args()

    df = pd.read_csv(args.path_to_run_report, header = None, names = ["name", "status", "exit", "workdir", "raw_GSE"])
    get_GSE_not_downloaded(df)

