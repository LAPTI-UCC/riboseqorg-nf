from numpy import inner
import xmltodict
import sys
import pandas as pd

# example = [{'@tag': 'cell line background', '#text': 'HCT-116'}, {'@tag': 'rna isolation', '#text': 'Total RNA'}, {'@tag': 'adapter sequence', '#text': 'TAGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'}]

def get_cell_info_from_dict_or_dict_list(list_or_dicts):
    """
    Given a dictionary or an unordered list of dictionaries (derived from the GEO report in xml format), returns the cell-line or strain of the sample.
    """
    cell_line = "No tag found"
    if type(list_or_dicts) == list:
        for dicts in list_or_dicts:
            if 'cell line' in dicts['@tag'] or "strain" in dicts['@tag'] or "tissue" in dicts['@tag']:
                cell_line = dicts['#text']
                return (cell_line)
            elif 'genotype' in dicts['@tag']:
                cell_line = dicts['#text']
                return (cell_line)
    elif type(list_or_dicts) == dict:         # code that checks the key
        if 'cell line' in list_or_dicts['@tag'] or "strain" in list_or_dicts['@tag'] or "tissue" in list_or_dicts['@tag']:
            cell_line = list_or_dicts['#text']
            return (cell_line)
        elif 'genotype' in dicts['@tag']:
            cell_line = dicts['#text']
            return (cell_line)

# 	GSE130465


def lister(Title, Organism, Cell_line, Description = "None_Available", Library_strategy = "None_Available", Protocol = "None_Available"):
    '''
    Just creates a list with the given arguments, allows for some of the to be optional.
    '''
    inner_list = [Title, Organism, Cell_line, Description, Library_strategy, Protocol]
    return (inner_list)


def parse_xml(xml_path):
    '''
    Reads the GSE*.xml into python and retrieves the fields needed, returned as a list within a dictionary, where each GSM is the key.
    '''
    output_dictionary = {}
    with open(xml_path) as xml_file:
        data_dict = xmltodict.parse(xml_file.read())

        # Progressively searches inside the MINiML, through each section, field, subfield, and sub_subfield for the data required.
        for MINiML in data_dict:
            for section in data_dict[MINiML]:
                if section == "Sample":
                    for field in data_dict[MINiML][section]:

                        # Set default values for those specifc sub-fields that are not always specified
                        desc = "No description available"
                        Lib_Strat = "Not available"
                        protocol = "Not available"

                        GSM_id = field["@iid"]

                        for subfield in field:
                            if subfield == "Title":
                                title = field[subfield]

                            if subfield == "Channel":
                                for sub_subfield in field[subfield]:
                                    
                                    # This snippet is left for a possibile future use.
                                    # if sub_subfield == "Source":
                                        # source_string = field[subfield][sub_subfield]

                                    if sub_subfield == "Extract-Protocol":
                                        protocol = field[subfield][sub_subfield]

                                    if sub_subfield == "Organism":
                                        organism = field[subfield][sub_subfield]['#text']

                                    if sub_subfield == "Characteristics":
                                        cell = get_cell_info_from_dict_or_dict_list(field[subfield][sub_subfield])
                            
                            if subfield == "Description":
                                desc = field[subfield]


                            if subfield == "Library-Strategy":
                                Lib_Strat = field[subfield]
                            
                        inner = lister(title,organism, cell, desc, Lib_Strat, protocol)
                        output_dictionary[GSM_id] = inner

    return ( output_dictionary )


def compile_df(dict):
    '''
    Creates a df from a dictionary. GSM are keys and relative fields are the values, given as a list.
    '''
    df = pd.DataFrame.from_dict(dict, orient="index")
    df.columns = ["Title", "Organism", "Cell/Strain", "Description", "Library_Strategy", "Extraction_Protocol"]
    return(df)

# Temporary functions, defined only for testing purposes and not required in the nextflow implementation.

def temp_read_path_list(txt_file):
    '''
    Reads a txt file containing the path to a GSE.xm. Each line of the file is a path and all paths are returned as a list.
    FUNCTION NEEDED ONLY FOR TESTING PURPOSES.
    '''
    file = open(txt_file, "r")
    path_list = file.readlines()
    return path_list

def temp_main(path_list):
    '''
    From a list of paths to GSE.xml reports, produces a csv file with all the single GSM and relative fields
    '''
    ls = [["","","","","",""]]
    output = pd.DataFrame(ls, columns = ["Title", "Organism", "Cell_line/Strain", "Description", "Library_Strategy", "Extraction_Protocol"])
    paths = temp_read_path_list(path_list)
    for GSE_report in paths:
        GSE_report = GSE_report[:-1]
        dict = parse_xml(GSE_report)
        df = compile_df(dict)
        output = pd.concat([output, df])
    output.to_csv("../../XML_CSV_output.csv")

# Function needed in the nextflow implementation, to name each signle csv accordingly.

def get_GSE(string):
    '''
    Using the input path, returns the GSE (which will be used to name the csv file appropriately).
    '''
    word_list = string.split("/")
    GSE = word_list[-1][:-4]
    return GSE


if __name__ == '__main__':
    #input_list = sys.argv[1]
    #temp_main(input_list)

    # When implemented in nextflow, this will be actual code:

    xml_path = sys.argv[1]
    dict = parse_xml(xml_path)
    df = compile_df(dict)
    GSE = get_GSE(xml_path)
    df.to_csv(GSE + ".csv")     # the final name should be something like "GSEnnnnnn_family.csv"


