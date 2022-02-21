import pandas as pd 



def convert_date_column(path_to_tsv, date_column):
    '''
    change date to yyyy/mm/dd from 'MONTH ABBR dd, yyyy'
    '''
    geo_ribosome_profiling = pd.read_csv(path_to_tsv)
    geo_ribosome_profiling.to_csv(path_to_tsv, sep="\t")

    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    for i in geo_ribosome_profiling.index:
        prev_date = geo_ribosome_profiling.at[i, f'{date_column}']
        month = prev_date.split(' ')[0]
        day = prev_date.split(' ')[1].strip(',')
        year = prev_date.split(' ')[2]
        geo_ribosome_profiling.at[i, f'{date_column}'] = year + "-" + str([i for i,x in enumerate(months) if x == month][0] + 1) + "-" + day

    geo_ribosome_profiling.to_csv(path_to_tsv, sep="\t")
