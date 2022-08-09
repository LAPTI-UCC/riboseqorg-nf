from numpy import inner
import xmltodict
import sys
import pandas as pd


def get_cell_info(tags_list):
    '''
    Given a dictionary or an unordered list of dictionaries (derived from the GEO report in xml format), returns the cell-line or strain of the sample.
    '''
    
    if type(tags_list) != list:
        tags_list = [tags_list]
    for dicts in tags_list:
        if 'cell' in dicts['@tag'] or "strain" in dicts['@tag']:
            return dicts['#text']
        elif "tissue" in dicts['@tag']:
            return dicts['#text']
        elif 'genotype' in dicts['@tag']:
            return dicts['#text']
    cell_line = "No tag found among 'Characteristics' detailing the cell/strain"
    return cell_line


def lister(Title, Organism, Source, Cell_line, Description = "None_Available", Library_strategy = "None_Available", Protocol = "None_Available", Tags = "None"):
    '''
    Creates a list with the given arguments, allows for some of the to be optional.
    '''
    inner_list = [Title, Organism, Source, Cell_line, Description, Library_strategy, Protocol, Tags]
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
                                    
                                    if sub_subfield == "Source":
                                        source_string = field[subfield][sub_subfield]

                                    if sub_subfield == "Extract-Protocol":
                                        protocol = field[subfield][sub_subfield]

                                    if sub_subfield == "Organism":
                                        organism = field[subfield][sub_subfield]['#text']

                                    if sub_subfield == "Characteristics":
                                        cell = get_cell_info(field[subfield][sub_subfield])
                                        
                                        # retrieves all the "characteristics tags" of the .xml 

                                        tags = field[subfield][sub_subfield]
                                        if type(tags) != list:
                                            tags = [tags]
                                        final_tags= {}
                                        for key in tags:
                                            new_key = key["@tag"]
                                            final_tags[new_key] = key["#text"]

                            if subfield == "Description":
                                desc = field[subfield]

                            if subfield == "Library-Strategy":
                                Lib_Strat = field[subfield]
                            
                        inner = lister(title,organism, source_string, cell, desc, Lib_Strat, protocol, final_tags)
                        output_dictionary[GSM_id] = inner

    return ( output_dictionary )


def compile_df(dict):
    '''
    Creates a df from a dictionary. GSM are keys and relative fields are the values, given as a list.
    '''
    df = pd.DataFrame.from_dict(dict, orient="index")
    df.columns = ["Title", "Organism", "Source", "Cell/Strain", "Description", "Library_Strategy", "Extraction_Protocol", "Tags"]
    return(df)


def get_GSE(string):
    '''
    Using the input path, returns the GSE (which will be used to name the csv file appropriately).
    '''
    word_list = string.split("/")
    GSE = word_list[-1][:-4]
    return GSE


if __name__ == '__main__':

    xml_path = sys.argv[1]
    dict = parse_xml(xml_path)
    df = compile_df(dict)
    GSE = get_GSE(xml_path)
    df.to_csv(GSE + ".csv")


