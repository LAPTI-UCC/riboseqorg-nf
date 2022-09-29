
import xmltodict
import sys
import pandas as pd


def get_sample_info(sample):
    '''
    Given a list of nested dictionaries, retrieves the required info and returns it in a tuple.
    '''

    # Set default values for returned argument in case they are missing

    GSM_id = sample["@iid"]
    desc = "No description available"
    Lib_Strat = "Not available"

    strain = ""
    cell = ""
    source = ""
    protocol = ""
    organism = ""
    final_tags = {}


    for section in sample:

        if section == "Title":
            title = sample["Title"]
        if section == "Description":
            desc = sample["Description"]
        if section == "Library-Strategy":
            Lib_Strat = sample["Library-Strategy"]
            
        if section == "Channel":
            for sub_section in sample["Channel"]:
                if sub_section == "Source":
                    source = sample["Channel"]["Source"]
                if sub_section == "Extract-Protocol":
                    protocol = sample["Channel"]["Extract-Protocol"]

# Improved: now this block of code accounts for multiple organisms in the same sample (ie. virus-infected human cells)
                if sub_section == "Organism":
                    organism_list = sample["Channel"]["Organism"]
                    for field in organism_list:
                        if "text" in field:
                            organism_list = [organism_list]
                    for org in organism_list:
                        organism = organism + " - " + org["#text"]
                    organism = organism[2:]

                if sub_section == "Characteristics":
                    cell = get_cell_info(sample[section][sub_section])
                    strain = get_strain_info(sample[section][sub_section])
                    tags = sample["Channel"]["Characteristics"]
                    if type(tags) != list:
                        tags = [tags]

                    for tag in tags:
                        if type(tag) == dict:
                            new_key = tag["@tag"]
                            final_tags[new_key] = tag["#text"]
                        elif type(tag) == str:
                            final_tags["tag1"] = tag 

    return (GSM_id, title, organism, source, strain, cell, desc, Lib_Strat, protocol, final_tags)


def get_cell_info(tags_list):
    '''
    Given a dictionary or an unordered list of dictionaries (derived from the GEO report in xml format), returns the cell line and tissue of the sample.
    '''
    if type(tags_list) != list:
        tags_list = [tags_list]
    cell_line = ""
    for tag in tags_list:
        get_tag_type = type(tag)
        if "str" not in str(get_tag_type):
            if 'cell' in tag['@tag']:
                cell_line = tag['#text']
            if 'tissue' in tag['@tag']:
                cell_line = cell_line + tag['#text']
        elif "str" in str(get_tag_type):
            cell_line = cell_line + " " + tag
    return cell_line


def get_strain_info(tags_list):
    '''
    Given a dictionary or an unordered list of dictionaries (derived from the GEO report in xml format), returns the strain and/or genotype of the sample.
    '''
    
    if type(tags_list) != list:
        tags_list = [tags_list]
    strain_info = ""
    for tag in tags_list:
        get_tag_type = type(tag)
        if "str" not in str(get_tag_type):
            if 'strain' in tag['@tag']:
                strain_info = tag['#text']
            if 'genotype' in tag['@tag']:
                strain_info = strain_info + tag['#text']
        elif "str" in str(get_tag_type):
            strain_info = strain_info + " " + tag
    return strain_info


def lister(Title, Organism, Source, Strain_Genotype, Cell_Tissue, Description, Library_strategy, Ribosome_position = "",Protocol = "None_Available", Tags = "None"):
    '''
    Creates a list with the given arguments, allows for some of the to be optional.
    '''
    inner_list = [Title, Organism, Source, Strain_Genotype, Cell_Tissue, Description, Library_strategy, Ribosome_position, Protocol, Tags]
    return inner_list


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

# Problem: if there is just one sample, then the nested structure is changed and the code will not work.
# Solution: I just need to put inside a list the study if it contains only one sample.
# Method: I check if data_dict[MINiML][section] contains traight up the ley "@iid" or not. If it does, then it's a single-sample study and needs to be put in
# a list. Otherwise, there are multiple samples and data_dict[MINiML][section] can be left as it is

                    # SHOULD I WRAP THIS UP IN A LIST?
                    sample_list = data_dict[MINiML][section]
                    for field in data_dict[MINiML][section]:
                        if field == "@iid":
                            sample_list = [data_dict[MINiML][section]]

                    for sample in sample_list:
                        GSM_id, title, organism, source, strain, cell, desc, Lib_Strat, protocol, final_tags = get_sample_info(sample)
                        rib_pos = ""
                        inner = lister(title,organism, source, strain, cell, desc, Lib_Strat, rib_pos, protocol, final_tags)
                        output_dictionary[GSM_id] = inner

    return  output_dictionary 


def compile_df(dictionary):
    '''
    Creates a df from a dictionary. GSM are keys and relative fields are the values, given as a list.
    '''
    #remove all new lines form the description section
    for key in dictionary:
        dictionary[key].append('')
        dictionary[key].append('')

        for idx, item in enumerate(dictionary[key]):
            if type(item) == str:
                dictionary[key][idx] = dictionary[key][idx].replace('\n', ' ').replace('\r', ' ')
    df = pd.DataFrame.from_dict(dictionary, orient="index")
    df.columns = ["Title", "Organism", "Source", "Strain/Genotype","Cell/Tissue", "Description", "Library_Strategy", "Ribosome_position", "Extraction_Protocol", "Tags", "Library_Strategy_Evidence", "Ribosome_Position_Evidence"]
    return df


def get_GSE(string):
    '''
    Using the input path, returns the GSE (which will be used to name the csv file appropriately).
    '''
    word_list = string.split("/")
    GSE = word_list[-1][:-4]
    return GSE


if __name__ == '__main__':

    xml_path = sys.argv[1]
    dictionary = parse_xml(xml_path)
    df = compile_df(dictionary)
    GSE = get_GSE(xml_path)
    df.to_csv(GSE + ".csv")


