import xmltodict
import sys


def parse_xml(xml_path):
    '''
    read xml into python 
    '''
    with open(xml_path) as xml_file:
        data_dict = xmltodict.parse(xml_file.read())
        for i in data_dict:
            for j in data_dict[i]:
                print(j)
                print(data_dict[i][j])
                print()


if __name__ == '__main__':
    xml_path = sys.argv[1]
    parse_xml(xml_path)
