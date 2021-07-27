import pandas as pd
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn

def get_raw_data(report_path):
    '''The function is used to generate raw data that is used in local data processing in servable module.'''

    smpl_info = pd.read_csv(report_path)
    smpl_info = smpl_info[smpl_info.sampling_date.str.match(r'[0-9]{4}-[0-9]{2}-[0-9]{2}')== True]
    
    unq_dict = {c_n:list(smpl_info[f'{c_n}'].drop_duplicates(keep='first')) for c_n in ['testing_lab','lineage', 'sampling_date']}
    unq_dict['sampling_date'] = sorted(unq_dict['sampling_date'])
    
    return smpl_info, unq_dict

def get_mut_data(report_path):
    '''The function is used to wrap metadata processing code.'''

    return pd.read_csv(report_path)

def quadratic_row_search(short_array, dataframe, variable, column_name):
    '''The function is used to generate a dataframe containing all rows from smpl_in that match query'''

    temp_df = pd.DataFrame()
    for record in short_array:
        query = (dataframe["sampling_date"] == str(record)) & (dataframe[str(column_name)] == str(variable))
        temp_df = temp_df.append(dataframe[query])

    return [variable, temp_df]


def quadratic_df_search(short_array, dataframe, variable, column_name, skip_null=False):
    '''The function is used to generates a list containing number of cases sorted by date matched to a variable in a returned list (from specified dataframe).'''
    temp_list = []

    for record in short_array:
        temp_list.append((record,len(dataframe[(dataframe["sampling_date"] == str(record)) & (dataframe[str(column_name)] == str(variable))])))
    
    temp_list.sort(key = lambda x:x[0])
    temp_list = [pair[1] for pair in temp_list]

    if skip_null: 
        if sum(temp_list)==0: temp_list=None

    return [variable, temp_list]









