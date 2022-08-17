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


def write_output_to_txt(GSEs_list, path = 'GSEs_not_downloaded.txt'):
    '''
    Given as input a list of GSEs, writes them to a txt file (whose path is given as second input).
    '''
    with open(path, 'w') as f:
        for GSE in GSEs_list:
            f.write(GSE)
            if GSEs_list.index(GSE) != len(GSEs_list)-1:
                f.write('\n')
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieves the GSEs that have not been downloaded correctly")
    parser.add_argument("path_to_report_run", type = str, help = "Path to the location of the run report, a csv file")
    parser.add_argument("txt_report_path", type = str, help = "Path the output txt file")
    args = parser.parse_args()

    df = pd.read_csv(args.path_to_report_run, header = None, names = ["name", "status", "exit", "workdir", "raw_GSE"])

    GSEs_list = get_GSE_not_downloaded(df)
    write_output_to_txt(GSEs_list, args.txt_report_path)