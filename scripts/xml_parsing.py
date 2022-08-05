from numpy import inner
import xmltodict
import sys
import pandas as pd

# example = [{'@tag': 'cell line background', '#text': 'HCT-116'}, {'@tag': 'rna isolation', '#text': 'Total RNA'}, {'@tag': 'adapter sequence', '#text': 'TAGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'}]

def get_cell_info_from_dict_list(list_of_dicts):
    """
    given an unordered list of dictionaries (derived from the GEO report in xml format), returns the cell-line or strain of the sample.
    """
    cell_line = "No tag found"
    def_type = "No tag found"
    for dicts in list_of_dicts:
        if 'cell line' in dicts['@tag'] or "strain" in dicts['@tag']:
            cell_line = dicts['#text']
            def_type = dicts['@tag']
    return (def_type, cell_line)


def lister(Title, Organism, Cell_line, Description = "None_Available", Library_strategy = "None_Available", Protocol = "None_Available"):
    inner_list = [Title, Organism, Cell_line, Description, Library_strategy, Protocol]
    return (inner_list)


def parse_xml(xml_path):
    '''
    read the GSE*.xml into python and retrieves the fields needed, returned as a dictionary for each GSM
    '''
    output_dictionary = {}
    with open(xml_path) as xml_file:
        data_dict = xmltodict.parse(xml_file.read())
        for i in data_dict:
            
            for j in data_dict[i]:
                if j == "Sample":
                    for k in data_dict[i][j]:
                        
                        # Set default values for those fields that are not always specified
                        desc = "No description available"
                        Lib_Strat = "Not available"
                        protocol = "Not available"

                        GSM_id = k["@iid"]

                        for each in k:
                            if each == "Title":
                                title = k[each]

                            if each == "Channel":
                                for every in k[each]:

                                    if every == "Source":

                                        source_string = k[each][every]

                                    
                                    if every == "Extract-Protocol":
                                        protocol = k[each][every]

                                    if every == "Organism":
                                        organism = k[each][every]['#text']

                                    if every == "Characteristics":
                                        s_type, cell = get_cell_info_from_dict_list(k[each][every])
                            
                            if each == "Description":
                                desc = k[each]

                            if each == "Library-Strategy":
                                Lib_Strat = k[each]
                            
                        inner = lister(title,organism, cell, desc, Lib_Strat, protocol)
                        output_dictionary[GSM_id] = inner
    return ( output_dictionary )


def compile_df(dict):
    '''
    Creates a df from a dictionary. GSM are keys and relative fields are the values, given as a list.
    '''
    df = pd.DataFrame.from_dict(dict, orient="index")
    df.columns = ["Title", "Organism", "Cell_line/Strain", "Description", "Library_Strategy", "Extraction_Protocol"]
    return(df)

def temp_read_path_list(txt_file):
    '''
    reads a txt file containing the path to a GSE.xm. Each line of the file is a path and all paths are returned as a list.
    FUNCTION NEEDED ONLY FOR TESTING PURPOSES.
    '''
    file = open(txt_file, "r")
    path_list = file.readlines()
    return path_list

def temp_main(path_list):
    '''
    from a list of paths to GSE.xml reports, produces a csv file with all the single GSM and relative fields
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


if __name__ == '__main__':
    input_list = sys.argv[1]
    temp_main(input_list)
    # output = 
    # for each in path_list:
        # 
        # 
        # output.append(df)
    # return output

    # xml_path = sys.argv[1]
    # dict = parse_xml(xml_path)
    # df = compile_df(dict)


