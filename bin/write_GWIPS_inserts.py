'''
This script takes in the path to the '*RunInfo.csv' file and the path to the
metadata csv file and writes the SQL inserts to add a study and all of its
tracks to GWIPS-viz.
The script assumes that the study metadata csv file has the following columns:
    - study_name
    - study_description
    - study_type
    - study_url
    - study_publication
This script also assumes that the inserts are for an existing Organism.
'''

import argparse
import pandas as pd
import sqlite3
import sys


def read_metadata(metadata_path):
    '''
    Reads the metadata csv file and returns a df of metadata.
    '''
    metadata = pd.read_csv(metadata_path)
    return metadata


def read_study_metadata(study_metadata_path):
    '''
    Reads the study metadata csv file and returns a df of study metadata.
    '''
    study_metadata = pd.read_csv(study_metadata_path)
    return study_metadata


def get_gwipsDB_for_organism(organism, annotationDB_path):
    '''
    Returns the path to the GWIPS-viz database for the given organism.
    '''
    conn = sqlite3.connect(annotationDB_path)
    c = conn.cursor()
    results = c.execute(
        f"SELECT * FROM gwips_organism WHERE scientific_name='{organism}';"
        ).fetchall()

    gwipsDB = results[0][1]
    return gwipsDB


def get_sample_prefix(study_metadata: pd.DataFrame) -> str:
    '''
    Generates the sample name prefix from the metadata df.

    returns a string of the sample name prefix
    '''
    first_author = study_metadata['authors'][0].split(',')[0].strip(' ')
    if len(first_author.split(" ")) > 1:
        first_author = first_author.split(" ")[1]

    date = study_metadata['date_published'][0].replace(' ', '_')
    sample_prefix = f"{first_author}_{date}"

    return sample_prefix


def get_library_type(metadata_row: pd.Series) -> str:
    '''
    Returns the library type based on the metadata row.
    '''
    if metadata_row[1]['Library_Strategy'] == 'RNA-Seq study':
        library_type = "mRNACov"

    elif metadata_row[1]['Library_Strategy'] == 'Ribo-seq study':
        if metadata_row[1]['Ribosome_position'] == 'Initiating':
            library_type = "RiboProInit"
        elif metadata_row[1]['Ribosome_position'] == 'Elongating':
            library_type = "RiboProElong"
        elif metadata_row[1]['Ribosome_position'] == 'Scanning':
            library_type = "RiboProScan"
        else:
            raise ValueError(f"Unknown ribosome position: {metadata_row[1]['Ribosome_position']}")
    else:
        raise ValueError(f"Unknown library strategy: {metadata_row['Library_Strategy']}")
    return library_type


def get_sample_long_label(metadata_row: pd.Series) -> str:
    '''
    Returns the sample long label based on the metadata row.
    '''
    library_type = get_library_type(metadata_row)

    long_label = f"Uniquely mapped reads from {metadata_row[1]['Cell/Tissue']} samples; {library_type};  {metadata_row[1]['Unnamed: 0']}"
    return long_label


def generate_sample_names(metadata: pd.DataFrame, study_metadata: pd.DataFrame) -> dict:
    '''
    Generates the sample names from the metadata df.

    returns a dictionary with GSE aas keys and sample name as the values
    '''
    sample_prefix = get_sample_prefix(study_metadata)

    sample_names = {}
    for row in metadata.iterrows():
        accession = row[1]['Unnamed: 0']
        library_suffix = get_library_type(row)
        if len(metadata) == 1:
            sample_name = f"{sample_prefix}_{library_suffix}"
        else:
            sample_name = f"{sample_prefix}{'_' + row[1]['Cell/Tissue']}_{library_suffix}"

        sample_names[accession] = sample_name
    return sample_names


def handle_sample_table(
        metadata: pd.DataFrame,
        sample_names: dict,
        gbdb_path: str
        ) -> list:
    '''
    Writes the SQL inserts to add the sample tables to GWIPS-viz.
    These tables store just the path to the sample bigWig file. The sample
    name is the table name and must match the trackDb entry.
    Returns a list of SQL inserts.
    '''

    sample_table_inserts = []

    for row in metadata.iterrows():
        accession = row[1]['Unnamed: 0']
        sample_name = sample_names[accession]
        sample_table_inserts.append(f"DROP TABLE IF EXISTS {sample_name};\n\n")
        sample_table_inserts.append(f"CREATE TABLE {sample_name} (filename VARCHAR(255));\n\n")
        sample_table_inserts.append(f"INSERT INTO {sample_name} VALUES ('{gbdb_path + sample_name}.bw');\n\n")

    return sample_table_inserts

def trackDb_parent_track_inserts(metadata: pd.DataFrame, sample_names: dict, study_metadata: pd.DataFrame) -> str:
    '''
    Writes the SQL insert to add the parent track to GWIPS-viz.
    Returns a string of the SQL insert.
    '''
    grp_dict = {
        'mRNACov': 'mRNA-Coverage',
        'RiboProInit': 'RP-InitiatingRibos',
        'RiboProElong': 'RP-ElongatingRibos',
        'RiboProScan': 'RP-Scanning',
    }

    color_dict = {
        'mRNACov': (0,0,200),
        'RiboProInit': (0,200,0),
        'RiboProElong': (200,0,0),
        'RiboProScan': (200,0,0),
    }
    parent_track_inserts = []

    # loop through metadata and make a parent track for each unique library type
    library_types = set()
    for row in metadata.iterrows():
        library_type = get_library_type(row)
        if library_type not in library_types:
            library_types.add(get_library_type(row))
            study_name = f"{get_sample_prefix(study_metadata)}_{grp_dict[library_type]}"

            parent_track_insert = f"INSERT INTO trackDb VALUES ('{study_name}', "
            parent_track_insert += f"'{get_sample_prefix(study_metadata).replace('_', ' ')}',"  # shortLabel
            parent_track_insert += f"'bigWig'," # type
            parent_track_insert += f"'{get_sample_long_label(row)}: {study_metadata['Title'][0]}',"  # longLabel
            parent_track_insert += f"0,"  # visibility
            parent_track_insert += f"1,"  # priority
            parent_track_insert += f"{color_dict[library_type][0]},"  # colorR
            parent_track_insert += f"{color_dict[library_type][1]},"  # colorG
            parent_track_insert += f"{color_dict[library_type][2]},"  # colorB
            parent_track_insert += f"{color_dict[library_type][0]},"  # altColorR
            parent_track_insert += f"{color_dict[library_type][1]},"  # altColorG
            parent_track_insert += f"{color_dict[library_type][2]},"  # altColorB
            parent_track_insert += "0,"  # useScore
            parent_track_insert += "0,"  # private
            parent_track_insert += "0,"  # restrictCount
            parent_track_insert += "'',"  # restrictList
            parent_track_insert += "'',"  # url
            html= f"""<h2>Description</h2> <p>{study_metadata['Title'][0]}.</p> <h2>Methods</h2>
                <p>Raw sequence data were obtained from NCBI GEO database (<a href="{study_metadata['GSE'][0]}">{study_metadata['Accession'][0]}</a>). Data from the following samples were processed:
                <table>
                <tr><th>Sample</th><th>Description</th></tr>
                """
            for row in metadata.iterrows():
                if get_library_type(row) == library_type:
                    html += f"<tr><td>{row[1]['Unnamed: 0']}</td><td>{row[1]['Title']}</td></tr>"
                else: 
                    print("FAIL")

            html += f"""
                </table></p>
                <p>Adaptor sequences had already been removed from reads. An alignment to ribosomal RNA was performed using <a href="http://bowtie-bio.sourceforge.net/index.shtml">Bowtie</a>, and aligning reads were discarded. An alignment to the mm10 genome assembly was then performed using <a href="http://bowtie-bio.sourceforge.net/index.shtml">Bowtie</a>, and these tracks contain the uniquely mapping reads that align to the A-site (using an offset of 15nt).
                </p>
                <h2>References</h2>
                <p>{study_metadata['authors'][0]} <a href="{study_metadata['doi'][0]}">.{study_metadata['title'][0]} </a>. <i>{study_metadata['journal'][0]}</i>.
                </p>
            """  # html

            parent_track_insert += f"'{html}',"  #html
            parent_track_insert += f"'{grp_dict[library_type]}'," #grp
            parent_track_insert += "1," #canPack
            parent_track_insert += "'compositeTrack on');\n\n" #settings
            parent_track_inserts.append(parent_track_insert)

    return parent_track_inserts


def write_trackDb_inserts(metadata: pd.DataFrame, sample_names: dict, study_metadata: pd.DataFrame) -> list:
    '''
    Writes the SQL inserts to add the trackDb entries to GWIPS-viz.
    Returns a list of SQL inserts.
    tableName | shortLabel  | type       | longLabel | visibility | priority | colorR | colorG | colorB | altColorR | altColorG | altColorB | 
    useScore | private | restrictCount | restrictList | url | html | grp      | canPack | settings
        '''

    grp_dict = {
        'mRNACov': 'mRNA-Coverage',
        'RiboProInit': 'RP-InitiatingRibos',
        'RiboProElong': 'RP-ElongatingRibos',
        'RiboProScan': 'RP-Scanning',
    }

    color_dict = {
        'mRNACov': (0,0,200),
        'RiboProInit': (0,200,0),
        'RiboProElong': (200,0,0),
        'RiboProScan': (200,0,0),
    }

    trackDb_inserts = []

    for row in metadata.iterrows():
        accession = row[1]['Unnamed: 0']
        sample_name = sample_names[accession]
        library_type = get_library_type(row)

        trackDb_insert = f"INSERT INTO trackDb VALUES ('{sample_name}'," #tableName
        trackDb_insert += f"'{get_sample_prefix(study_metadata).replace('_', ' ')}'," #shortLabel
        trackDb_insert += f"'bigWig'," #type
        trackDb_insert += f"'{get_sample_long_label(row)}'," #longLabel
        trackDb_insert += f"1," #visibility
        trackDb_insert += f"1," #priority
        trackDb_insert += f"{color_dict[library_type][0]}," #colorR
        trackDb_insert += f"{color_dict[library_type][1]}," #colorG
        trackDb_insert += f"{color_dict[library_type][2]}," #colorB
        trackDb_insert += f"{color_dict[library_type][0]}," #altColorR
        trackDb_insert += f"{color_dict[library_type][1]}," #altColorG
        trackDb_insert += f"{color_dict[library_type][2]}," #altColorB
        trackDb_insert += f"0," #useScore
        trackDb_insert += f"0," #private
        trackDb_insert += f"0," #restrictCount
        trackDb_insert += f"''," #restrictList
        trackDb_insert += f"''," #url
        trackDb_insert += f"''," #html
        trackDb_insert += f"'{grp_dict[library_type]}'," #grp
        trackDb_insert += f"1," #canPack
        trackDb_insert += f"'autoScale on\nalwaysZero on\nparent {get_sample_prefix(study_metadata)}_{grp_dict[library_type]} off');" #settings
        trackDb_inserts.append(trackDb_insert)

    return trackDb_inserts



def write_study_inserts(metadata: pd.DataFrame, sample_names: dict, gbdb_path: str="/gbdb/") -> list:
    '''
    Writes the SQL inserts to add a study to GWIPS-viz.
    '''
    study_inserts = []

    study_inserts = handle_sample_table(metadata, sample_names, gbdb_path)

    return study_inserts


def subset_metadata(metadata: pd.DataFrame, study_metadata: pd.DataFrame) -> pd.DataFrame:
    '''
    Subset metadata to only include samples that have non NA values for the following columns:
    - Organism
    - Library type
    - Ribosome Position
    Returns a subsetted metadata dataframe.
    '''
    subset = metadata[metadata['Organism'].notna()]
    subset = subset[subset['Library_Strategy'].notna()]
    subset = subset[subset['Ribosome_position'].notna()]
    return subset


def main(args):
    '''
    Main function.
    '''
    metadata = read_metadata(args.m)
    study_metadata = read_study_metadata(args.s)
    metadata = subset_metadata(metadata, study_metadata)

    if metadata.empty:
        print("No valid samples found in metadata. Valid entries must be present for 'Organism', 'Library_Strategy' and 'Ribosome_position'.Exiting.")
        sys.exit(1)


    #loop through unique organisms in study_metadata and create output insert files for each
    for organism in metadata['Organism'].unique():
        organism_metadata = metadata[metadata['Organism'] == organism]

        organism = organism.strip(' ')

        with open(f"{get_gwipsDB_for_organism(organism, args.db)}.sql", "w") as f:
            sample_names = generate_sample_names(organism_metadata, study_metadata)

            # Handle sample tables with paths to bigWig files
            study_inserts = write_study_inserts(organism_metadata, sample_names)

            f.write(f"--\n -- {get_sample_prefix(study_metadata)}\n--\n")
            for study_insert in study_inserts:
                f.write(study_insert)

            # Handle trackDb inserts. Each sequence type in each study gets  a parent track and a track for each sample
            f.write("--\n -- trackDb Parent Tracks\n--\n")
            parent_tracks = trackDb_parent_track_inserts(organism_metadata, sample_names, study_metadata)

            for parent_track in parent_tracks:
                f.write(parent_track)

            f.write("--\n -- trackDb Child Tracks\n--\n")
            trackDb = write_trackDb_inserts(organism_metadata, sample_names, study_metadata)
            for track in trackDb:
                f.write(track)

    return study_inserts



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Writes the SQL inserts to add a study and all of its tracks to GWIPS-viz")
    parser.add_argument("-s", type = str, help = "Path to the location of the study metadata csv file")
    parser.add_argument("-m", type = str, help = "Path to the location of the metadata csv file")
    parser.add_argument("-g", type = str, help = "Path to the base gbdb directory", default="/gbdb/")
    parser.add_argument("--db", type = str, help = "Path to the annotation inventory sqlite database")
    args = parser.parse_args()
    main(args)


'''
python scripts/write_GWIPS_inserts.py -s /home/jack/projects/Riboseq-Database/data/GSE73136_metadata.csv -m /home/jack/projects/Riboseq-Database/null/null/GSE73136_family.csv --db /home/jack/projects/riboseq_data_processing/annotation_inventory/annotation_inventory.sqlite
python scripts/write_GWIPS_inserts.py -s /home/jack/projects/Riboseq-Database/data/GSE115161_metadata.csv -m /home/jack/projects/Riboseq-Database/null/null/GSE115161_family.csv --db /home/jack/projects/riboseq_data_processing/annotation_inventory/annotation_inventory.sqlite
'''