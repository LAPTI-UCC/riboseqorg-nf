import argparse
import numpy as np
import pandas as pd


def load_csv(path):
    '''
    Loads the .csv, given as an argument to the script in the terminal, as a dataframe (df).
    '''

    df = pd.read_csv(path, )
    
    return(df)


def count_empty_cells_library_strategy(df):
    '''
    Check a provided column and count the number of empty cells
    '''
    count = 0
    for idx, row in df.iterrows():
        if row['Library_Strategy'] in ['Ribo-seq study', 'RNA-seq study']:
            count += 1
    num_rows = len(df.index)
    prop = count/num_rows
    if prop != 1:
        return Exception (f"Incomplete Library Strategy. {round(prop *100,2)}% complete")
    
    return prop


def main(args):
    '''
    Run the column checker and make decisions on completeness. Raise errors where needed
    '''

    df = load_csv(args.Path)
    ls_check = count_empty_cells_library_strategy(df)
    if type(ls_check) == Exception:
        print(args.Path)
        raise ls_check


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Checks how well metadata collection went")
    parser.add_argument("Path", type = str, help = "Path to metadata report csv")
    args = parser.parse_args()

    main(args)