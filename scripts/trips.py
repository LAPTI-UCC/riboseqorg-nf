"""
This script is used to add a study to the trips.sqlite database. 

The script requires access to the trips.sqlite database. The database is located on the prod server.
The insert statements can be written to a file and then run on the prod server but this requires a 
clone of the DB so that the id numbers are consistent. 

The other option is to use port forwarding to access the database directly. This is probably the easiest way to do it. But I haven't worked it out yet.

Alternatively, the script can be run on the prod server and interact with the database directly. This is the simplest but probably not the best way to do it.

Arguments:
    -s, --study: Path to study metadata file
    -r, --run: Path to SRA runInfo file
    -m, --mode: 'db' or 'file'. If 'file', the script will write to a file. If 'db', the script will write to the database directly.
    -d, --database: Path to trips.sqlite database (or clone if in file mode)
    -o, --output: Path to output file for file mode

"""

import sqlite3
import argparse
import pandas as pd


def get_study_metadata(study_metadata_file: str) -> dict:
    """
    Return a dictionary of study metadata

    Arguments:
        study_metadata_file: <str> Path to the study metadata file (CSV)

    Returns:
        study_metadata: <dict> Dictionary of study metadata
    """
    study_metadata = pd.read_csv(
        study_metadata_file, sep=",", header=0, index_col=False
    )
    study_metadata_dict = {
        key: study_metadata[key][0] for key in study_metadata.to_dict()
    }

    return study_metadata_dict

def get_run_df(run_info_file: str) -> pd.DataFrame:
    """
    Return a pandas dataframe of the run info file

    Arguments:
        run_info_file: <str> Path to the SRA runInfo file

    Returns:
        run_df: <pd.DataFrame> pandas dataframe of the run info file
    """
    run_df = pd.read_csv(
        run_info_file,
        header=0,
        index_col=False,
    )

    return run_df


def get_new_study_id(db: str) -> int:
    """
    Return a new study id to enable inserts into the db

    Arguments:
        db: <str> Path to the sqlite database

    Returns:
        new_key: <int> new integer key for the study that is not used in the studies database
    """
    connection = sqlite3.connect(db)
    cursor = connection.cursor()
    result = cursor.execute("SELECT MAX(study_id) FROM studies;").fetchall()[0][0]
    new_key = int(result) + 1

    return new_key


def get_new_file_id(db: str) -> int:
    """
    Return a new file id to enable inserts into the db

    Arguments:
        db: <str> Path to the sqlite database

    Returns:
        new_key: <int> new integer key for the file that is not used in the files database
    """
    connection = sqlite3.connect(db)
    cursor = connection.cursor()
    result = cursor.execute("SELECT MAX(file_id) FROM files;").fetchall()[0][0]
    new_key = int(result) + 1

    return new_key


def get_organism_id(organism: str, db: str) -> int:
    """
    Return the organism id for the organism of interest.

    Arguments:
        organism: <str> Name of the organism
        db: <str> Path to the sqlite database

    Returns:
        organism_id: <int> organism id for the organism of interest
    """
    connection = sqlite3.connect(db)
    cursor = connection.cursor()
    result = cursor.execute(
        f"SELECT organism_id FROM organisms WHERE organism_name = '{organism}';"
    ).fetchall()[0][0]
    organism_id = int(result)

    return organism_id


def get_organism_studyID_dict(study_metadata: pd.DataFrame, db: str) -> dict:
    """
    Return a dictionary of organism names and study ids

    Arguments:
        study_metadata: <pd.DataFrame> pandas dataframe of the study metadata
        db: <str> Path to the sqlite database

    Returns:
        organism_studyID_dict: <dict> Dictionary of organism names and study ids
    """
    organism_studyID_dict = {}
    study_id = get_new_study_id(db)

    for organism in study_metadata["Organism"].split(";"):
        organism = organism.strip().replace(" ", "_").lower()
        organism_id = get_organism_id(organism, db)
        organism_studyID_dict[organism_id] = study_id
        study_id += 1

    return organism_studyID_dict

def handle_study(study_metadata: dict, db: str, study_organism_dict: dict) -> list:
    """
    Handle the inserts into the study table of the trips data base.

    table nam is `studies` with the following fields:
        `study_id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
        `organism_id` integer DEFAULT NULL
        `study_name` varchar(500) DEFAULT NULL
        `paper_authors` varchar(300) DEFAULT NULL
        `srp_nos` varchar(300) DEFAULT NULL
        `paper_year` varchar(300) DEFAULT NULL
        `paper_pmid` varchar(300) DEFAULT NULL
        `paper_link` varchar(300) DEFAULT NULL
        `gse_nos` varchar(300) DEFAULT NULL
        `adapters` varchar(300) DEFAULT NULL
        `paper_title` varchar(300) DEFAULT NULL
        `description` varchar(2000) DEFAULT NULL
        `private` integer DEFAULT NULL
        `owner` integer DEFAULT NULL

    Parameters:
        study_metadata: <dict> Dictionary of study metadata
        db: <str> Path to the sqlite database
        study_organism_dict: <dict> Dictionary of organism names and study ids

    Returns:
        statements: <list> List of insert statements for the study table

    """
    statements = []

    for organism in study_metadata["Organism"].split(";"):
        organism = organism.strip().replace(" ", "_").lower()
        organism_id = get_organism_id(organism, db)
        title = study_metadata["Title"]
        authors = study_metadata["authors"]
        srp_nos = study_metadata["SRA"]
        paper_year = study_metadata["date_published"]
        paper_pmid = study_metadata["PMID"]
        paper_link = f"http://doi.org/{study_metadata['doi']}"
        gse_nos = study_metadata["Accession"]
        adapters = ""  # adatpers are not included in the study metadata file so this is left blank for now until I incorporate them
        paper_title = study_metadata["title"]
        description:str = study_metadata["abstract"]
        private = 0
        owner = 1

        insert_statement = f"INSERT INTO studies VALUES ({study_organism_dict[organism_id]}, {organism_id}, '{title}', '{authors}', '{srp_nos}', '{paper_year}', '{paper_pmid}', '{paper_link}', '{gse_nos}', '{adapters}', '{paper_title}', '{description}', {private}, {owner});"
        statements.append(insert_statement)

    return statements


def handle_files(run_info_file: str, db: str, study_organism_dict: dict) -> list:
    """
    Handle the inserts into the files table of the trips data base.

    Parameters:
        run_info_file: <str> Path to the SRA runInfo file
        db: <str> Path to the sqlite database
        study_organism_dict: <dict> Dictionary of organism ids (keys) and study ids (values)
    
    Returns:
        statements: <list> List of insert statements for the files table

    `files` (
  `file_id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `organism_id` integer DEFAULT NULL
,  `study_id` integer DEFAULT NULL
,  `file_name` varchar(300) DEFAULT NULL
,  `file_description` varchar(300) DEFAULT NULL
,  `file_type` varchar(300) DEFAULT NULL
,  `owner` integer DEFAULT NULL
, mapped_reads INT, control TINYINT, cell_line VARCHAR);
    """
    statements = []
    run_df = get_run_df(run_info_file)

    file_id = get_new_file_id(db)

    for organism in run_df["Organism"].unique():
        organism_df = run_df[run_df["Organism"] == organism]
        organism = organism.strip().replace(" ", "_").lower()
        study_id = study_organism_dict[get_organism_id(organism, db)]
        for index, row in organism_df.iterrows():
            file_name = row["Title"]
            file_description = row["Source"]
            if type(row['Library_Strategy']) != str:
                continue
                # raise ValueError("Library_Strategy is not defined for this file")
            if "ribo" in row['Library_Strategy'].lower():
                file_type = "riboseq"
            elif "rna" in row['Library_Strategy'].lower():
                file_type = "rnaseq"
            else:
                file_type = "unknown"
            owner = 1
            mapped_reads = 0
            control = 0
            cell_line = ""

            insert_statement = f"INSERT INTO files VALUES ({file_id}, {get_organism_id(organism, db)}, {study_id}, '{file_name}', '{file_description}', '{file_type}', {owner}, {mapped_reads}, {control}, '{cell_line}');"
            statements.append(insert_statement)
            file_id += 1

    return statements


def main(args):
    """
    Main function that runs the script depending on the mode specified
    """
    study_metadata = get_study_metadata(args.study)

    study_organism_dict = get_organism_studyID_dict(study_metadata, args.database)
    study_statements = handle_study(study_metadata, args.database, study_organism_dict)

    file_statements = handle_files(args.run, args.database, study_organism_dict)
    statements = study_statements + file_statements
    if args.mode == "file":
        with open(args.output, "w") as f:
            for statement in statements:
                f.write(statement)

    elif args.mode == "db":
        connection = sqlite3.connect(args.database)
        cursor = connection.cursor()
        for statement in statements:
            cursor.execute(statement)
        connection.commit()
        connection.close()
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--study", help="Path to study metadata file", required=True
    )
    parser.add_argument("-r", "--run", help="Path to SRA runInfo file", required=True)
    parser.add_argument(
        "-m", "--mode", help="db or file", choices=["db", "file"], required=True
    )
    parser.add_argument(
        "-d",
        "--database",
        help="Path to trips.sqlite database (or clone if in file mode)",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", help="Path to output file for file mode", required=False
    )
    args = parser.parse_args()
    main(args)


"""
 python scripts/write_GWIPS_inserts.py -s /home/jack/projects/Riboseq-Database/data/GSE115161_metadata.csv 
                                        -m /home/jack/projects/Riboseq-Database/null/null/GSE115161_family.csv 
                                        --db /home/jack/projects/riboseq_data_processing/annotation_inventory/annotation_inventory.sqlite


python scripts/trips.py -s /home/jack/projects/Riboseq-Database/data/GSE115161_metadata.csv -r /home/jack/projects/Riboseq-Database/null/null/GSE115161_family.csv -m db -d /home/jack/projects/Riboseq-Database/trips.sqlite


CREATE TABLE `files` (
  `file_id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `organism_id` integer DEFAULT NULL
,  `study_id` integer DEFAULT NULL
,  `file_name` varchar(300) DEFAULT NULL
,  `file_description` varchar(300) DEFAULT NULL
,  `file_type` varchar(300) DEFAULT NULL
,  `owner` integer DEFAULT NULL
, mapped_reads INT, control TINYINT, cell_line VARCHAR);

CREATE TABLE sqlite_sequence(name,seq);

CREATE TABLE `organisms` (
  `organism_id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `organism_name` varchar(300) DEFAULT NULL
,  `transcriptome_list` varchar(500) DEFAULT NULL
,  `gwips_databasename` varchar(300) DEFAULT NULL
,  `gwips_clade` varchar(300) DEFAULT NULL
,  `gwips_organism` varchar(300) DEFAULT NULL
,  `gwips_database` varchar(300) DEFAULT NULL
,  `default_transcript` varchar(300) DEFAULT NULL
,  `private` integer DEFAULT NULL
, owner INT(1));

CREATE TABLE `studies` (
  `study_id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `organism_id` integer DEFAULT NULL
,  `study_name` varchar(500) DEFAULT NULL
,  `paper_authors` varchar(300) DEFAULT NULL
,  `srp_nos` varchar(300) DEFAULT NULL
,  `paper_year` varchar(300) DEFAULT NULL
,  `paper_pmid` varchar(300) DEFAULT NULL
,  `paper_link` varchar(300) DEFAULT NULL
,  `gse_nos` varchar(300) DEFAULT NULL
,  `adapters` varchar(300) DEFAULT NULL
,  `paper_title` varchar(300) DEFAULT NULL
,  `description` varchar(2000) DEFAULT NULL
,  `private` integer DEFAULT NULL
,  `owner` integer DEFAULT NULL
);

"""
