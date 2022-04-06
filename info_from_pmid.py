import subprocess


def multiline_parse(medline, idx):
    '''
    check if following lines are of the same entry. Merge enteries that span multi lines and returna as string.
    '''
    return_string = medline[idx].split('- ')[1]
    
    cont = True
    while cont:
        if medline[idx] == medline[-1]:
            cont = False
        else:
            if len(medline[idx + 1].split("- ")) == 2:
                cont = False
                continue 
            else:
                return_string += medline[idx + 1] 
                idx += 1
    
    return ' '.join(return_string.split())


def get_info_from_medline(query):
    '''
    take a query (doi or pmid work) and return paper metadata in medline format and further parse this into an information dictionary 
    '''
    info_dict = {
        "PMID":None,
        "authors":None,
        "abstract":None,
        "title":None,
        "doi":None,
        "date_published":None,
        "PMC":None,
        "journal": None
    }
    if ";" in str(query):
        query = query.split(';')[0]

    s = subprocess.check_output(f'esearch -db pubmed -query {query} | efetch -format medline', stderr=subprocess.STDOUT, shell=True)
    medline = str(s).split('\\n')
    
    for i in enumerate(medline):
        line = i[1].split("- ")
        if len(line) == 2:
            entry = multiline_parse(medline, i[0])

            if "AB" in line[0]: 
                info_dict['abstract'] = entry 

            if "PMID" in line[0]: 
                info_dict['PMID'] = entry 

            if "JT" in line[0]: 
                info_dict['journal'] = entry

            if "TI" in line[0]: 
                info_dict['title'] = entry 

            if "PMC" in line[0]: 
                info_dict['PMC'] = entry

            if "DP" in line[0]: 
                info_dict['date_published'] = entry

            if "AID" in line[0] and "doi" in line[1]: 
                info_dict['doi'] = entry.split(' ')[0]
                
            if "FAU" in line[0] and info_dict['authors'] == None:
                name = entry.split(',')
                info_dict['authors'] = ' '.join(name[1:]) + " " + name[0] + ", "
            elif "FAU" in line[0] and info_dict['authors'] != None:
                name = entry.split(',')
                info_dict['authors'] += ' '.join(name[1:]) + " " + name[0] + ", "
    try:        
        info_dict['authors'] = info_dict['authors'][:-2]
    except:
        info_dict['authors'] = info_dict['authors']
    return info_dict
