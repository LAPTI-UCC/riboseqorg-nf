import argparse
import pandas as pd



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Writes the SQL inserts to add a study and all of its tracks to GWIPS-viz")
    parser.add_argument("-s", type = str, help = "Path to the location of the '*RunInfo.csv' file")
    parser.add_argument("-s", type = str, help = "Path to the location of the metadata csv file")

    args = parser.parse_args()

