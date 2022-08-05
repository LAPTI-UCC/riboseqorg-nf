import xmltodict
import sys

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

def inner_dictionary_compiler(Title, Organism, Cell_line, Description, Library_strategy, Protocol = ""):
    '''
    Compiles a dictionary with the relevant fields for the processed GSM.
    '''
    inner_dic = {"Title": Title, "Organism": Organism, "Cell_line/Strain" : Cell_line, "Description" : Description, "Library_strategy" : Library_strategy, "Extraction_protocol" : Protocol}
    return inner_dic

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
                        #print()
                        #print("-"*40)
                        GSM_id = k["@iid"]
                        #print("ID: ", GSM_id)
                        #print("-" *40)
                        for each in k:
                            if each == "Title":
                                title = k[each]
                                #print("title: ", title)
                            if each == "Channel":
                                for every in k[each]:

                                    if every == "Source":
                                        #print (every)
                                        #print( "~" * 50)
                                        source_string = k[each][every]
                                        #print(source_string)
                                    
                                    if every == "Extract-Protocol":
                                        protocol = k[each][every]
                                        #print (every)
                                        #print( "~" * 50)
                                        #print(protocol)
                                    if every == "Organism":
                                        organism = k[each][every]['#text']
                                        #print (every)
                                        #print( "~" * 50)
                                        #print(organism)

                                    if every == "Characteristics":
                                        s_type, cell = get_cell_info_from_dict_list(k[each][every])
                                        #print (every)
                                        #print( "~" * 50)
                                        #print(s_type, ": ", cell)
                                    #print()
                                #print()
                            #print()
                            
                            if each == "Description":
                                desc = k[each]
                                #print(each)
                                #print("~"*50)
                                #print(desc)
                                #print()
                            if each == "Library-Strategy":
                                Lib_Strat = k[each]
                                #print(each)
                                #print("~"*50)
                                #print(Lib_Strat)

                                        #cell_line = k[each][every]['#text']
                                        #strategy = k[each][every]['#text']
                                    #print(k[each][every])

                        inner_dic = inner_dictionary_compiler(title,organism, cell, desc, Lib_Strat)
                        output_dictionary[GSM_id] = inner_dic
    return output_dictionary





if __name__ == '__main__':
    xml_path = sys.argv[1]
    parse_xml(xml_path)


