#Standard libraries
import concurrent.futures, pandas as pd, os, requests, csv, subprocess, time, sys
from pathlib import Path
from os.path import dirname, join
from datetime import date, datetime, timedelta
#Bokeh modules
from bokeh.io.doc import curdoc
from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource
from bokeh.models.layouts import Column
from bokeh.models.widgets import Button, Select, DateRangeSlider, DataTable, TableColumn 
from bokeh.models.widgets.groups import CheckboxGroup
#Own modules
from get_data import get_raw_data, get_mut_data, quadratic_row_search, quadratic_df_search
from plotting_functions import stacked_bar_plot



#General data
report_time = time.strftime("%m_%d_%Y")
report_path = sys.argv[1]


full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[1]
raw_data = get_raw_data(report_path)
mut_data = get_mut_data(report_path)


['receiving_lab_sample_id','testing_lab','testing_lab_sample_id','sample_type','normalized_sample_type','sampling_date','notification_date','seq_date','seq_institution',
'result','ordering_institution','E_gene','ORF1a','RdRP/S','lineage','genome_length','genome_N_percentage','genome_GC_content','SEQ_DATE',
'COVERAGE(%)','AVERAGE_COVERAGE','MAPPED_FRACTION','READS_MAPPED','TOTAL_READS']
lv_name_list = [
    'Parauga ID',
    'Testēšanas laboratorija',
    'Par. veids',
    'Par. ņemšanas datums']

metadata_list = [
    'receiving_lab_sample_id',
    'testing_lab',
    'normalized_sample_type',
    'sampling_date']

# lv_name_list = [
#     'Parauga ID',
#     'Testēšanas laboratorija',
#     'Par. veids',
#     'Par. ņemšanas datums',
#     'Rezultātu paziņošanas datums',
#     'Sekvenēšanas datums',
#     'Sekvenēšanas vieta',
#     'RLPĶR Rezultāts',
#     'Sekvencēšanas pasūtītājs',
#     'Ct - E-gēns',
#     'Ct - ORF1a',
#     'Ct - RdRP/S',
#     'Celms',
#     'Consensus sekv. garums',
#     'Consensus sekv. N saturs',
#     'Consensus sekv. GC saturs',
#     'Genoma pārklājums(%)',
#     'Vidējais pārklājums',
#     'Kartēta lasījumu daļa(%)',
#     'Kartēto lasījumu skaits',
#     'Kopējais lasījumu skaits']

#Data table
# metadata_list = [
#     'receiving_lab_sample_id',
#     'testing_lab',
#     'normalized_sample_type',
#     'sampling_date',
#     'notification_date',
#     'seq_date',
#     'seq_institution',
#     'result',
#     'ordering_institution',
#     'E_gene',
#     'ORF1a',
#     'RdRP/S',
#     'lineage',
#     'genome_length',
#     'genome_N_percentage',
#     'genome_GC_content',
#     'COVERAGE(%)',
#     'AVERAGE_COVERAGE'
#     ,'MAPPED_FRACTION'
#     ,'READS_MAPPED'
#     ,'TOTAL_READS'
# ] 
mutation_list = ['B.1.1.7+E484K','B.1.616','B.1.617','CLUSTER_5','E484K','N501Y','ORF1a(del3675-3677)','P.1','S_GENE_DELETION','Y453F','B.1.525','B.1.427/B.1.429','P.3','B.1.617.1','B.1.617.2','B.1.617.3','B.1.620','B.1.621','B.1.351',"B.1.1.7"]
combined_data = mut_data #for inhouse data
Columns = [TableColumn(field=Ci, title=Ci) for Ci in lv_name_list+mutation_list] # bokeh columns
rename_dict = {eng:lv for eng,lv in zip(metadata_list,lv_name_list)}

#Extracting sample_by_date data
smpl_info, unq_dict = raw_data[0], raw_data[1]
smpl_by_date = [len(smpl_info.query(f"sampling_date == '{date}'")) for date in unq_dict['sampling_date']]

#Extracting dataframes that refer to a single testing_lab
by_ldf_dict = {lab:pd.DataFrame() for lab in unq_dict['testing_lab']}
with concurrent.futures.ProcessPoolExecutor() as executor: # applying the function in-parallel on different cores
    results = [executor.submit(quadratic_row_search, unq_dict['sampling_date'], smpl_info, lab, 'testing_lab') for lab in unq_dict['testing_lab']]#submitting function calls to different processes
    for f in concurrent.futures.as_completed(results): #collecting processing result
        by_ldf_dict[f.result()[0]] = f.result()[1]

#For each testing lab, extracting lineage distribution for each sampling date
by_lab_lin_dict = {}
for key in by_ldf_dict.keys():
    by_lab_lin_dict[key] = {}
    with concurrent.futures.ProcessPoolExecutor() as executor: # applying the function in-parallel on different cores
        results = [executor.submit(quadratic_df_search, unq_dict['sampling_date'], by_ldf_dict[key], lin, 'lineage', True) for lin in unq_dict['lineage']]#submitting function calls to different processes
        for f in concurrent.futures.as_completed(results): #collecting processing result
            if f.result()[1] is not None: by_lab_lin_dict[key][f.result()[0]] = f.result()[1]
            else:continue
                
#Generating sample_by_date_by_lab_by_lin relative data
by_lab_lin_dict_rel = {key:{lin:[round(100*i/j,2) for i,j in zip(by_lab_lin_dict[key][lin], smpl_by_date)] for lin in by_lab_lin_dict[key].keys()} for key in by_lab_lin_dict.keys()}
for key in by_lab_lin_dict.keys():
    by_lab_lin_dict[key]['sampling_date'] = unq_dict['sampling_date'] #adding sampling data to the data dict
    by_lab_lin_dict[key]['total_cases'] = smpl_by_date #add total number of samples to the data dict
    by_lab_lin_dict_rel[key]['sampling_date'] = unq_dict['sampling_date'] #adding sampling data to the data dict
    by_lab_lin_dict_rel[key]['total_cases'] = smpl_by_date #add total number of samples to the data dict

by_lin_dict = {}
with concurrent.futures.ProcessPoolExecutor() as executor: # applying the function in-parallel on different cores
    results = [executor.submit(quadratic_df_search, unq_dict['sampling_date'], smpl_info, lin, 'lineage') for lin in unq_dict['lineage']]#submitting function calls to different processes
    for f in concurrent.futures.as_completed(results): #collecting processing results to the dict
        by_lin_dict[f.result()[0]] = f.result()[1]
by_lin_dict_rel = {key:[round(100*i/j,2) for i,j in zip(by_lin_dict[key], smpl_by_date)] for key in by_lin_dict.keys()}
by_lin_dict_rel['sampling_date'] = unq_dict['sampling_date'] #adding sampling data to the data dict
by_lin_dict_rel['total_cases'] = smpl_by_date #add total number of samples to the data dict

p_abs_list = {lab:stacked_bar_plot(by_lab_lin_dict[lab],unq_dict,'Celms') for lab in by_lab_lin_dict.keys()}
by_lab_lin_dict['Kopā'] = by_lin_dict
by_lab_lin_dict['Kopā']['sampling_date'] = unq_dict['sampling_date'] #adding sampling data to the data dict
by_lab_lin_dict['Kopā']['total_cases'] = smpl_by_date
p_abs_list['Kopā'] = stacked_bar_plot(by_lin_dict,unq_dict,'Celms')

p_rel_list = {lab:stacked_bar_plot(by_lab_lin_dict_rel[lab],unq_dict,'Celms', units='rel') for lab in by_lab_lin_dict_rel.keys()}
by_lab_lin_dict_rel['Kopā'] = by_lin_dict_rel
by_lab_lin_dict_rel['Kopā']['sampling_date'] = unq_dict['sampling_date'] #adding sampling data to the data dict
by_lab_lin_dict_rel['Kopā']['total_cases'] = smpl_by_date
p_rel_list['Kopā'] = stacked_bar_plot(by_lin_dict_rel,unq_dict,'Celms', units = 'rel')

#Data wrangling to download
data_to_download = ColumnDataSource(data=smpl_info)

# #Creating bokeh document to serve
bokeh_doc = curdoc()

#Custom interactivity
global lg_r_list #declaring global variable to store legend renderers as they get updated
def legend_callback(attr, old, new):
    range_start = datetime.fromtimestamp(r_slider.value[0] / 1e3) #get datetime object from slider current position - range start
    range_end = datetime.fromtimestamp(r_slider.value[1] / 1e3) #get datetime object from slider current position - range end
    list_of_dates = date_range(range_start,range_end)
    list_of_visible_lg = [lg for lg in legend if curdoc().get_model_by_name(lg) and curdoc().get_model_by_name(lg).visible]
    if curdoc().get_model_by_name('select').value == 'Kopā':
        df = combined_data.loc[(combined_data['sampling_date'].isin(list_of_dates)) & (combined_data['lineage'].isin(list_of_visible_lg))]
        df = df.rename(columns = rename_dict)
    else:
        df = combined_data.loc[(combined_data['testing_lab'] == curdoc().get_model_by_name('select').value) & (combined_data['sampling_date'].isin(list_of_dates)) & (combined_data['lineage'].isin(list_of_visible_lg))]
        df = df.rename(columns = rename_dict) #rename columns to lv
    data_table = DataTable(columns=Columns, source=ColumnDataSource(df), autosize_mode="fit_columns", aspect_ratio = 1, align = 'end', name = 'table') #generate bokeh table
    curdoc().get_model_by_name('row').children[1].children[0] = data_table #update bokeh table in the layout
    
labs = [key for key in p_abs_list.keys()] #used to select lab from drop-down list
empty_fig = figure(plot_width = 1000, plot_height = 1000, name='empty') #placeholder when no figure is selected

select = Select(title = "Testēšanas laboratorija:", value = "", options = [""] + labs, name = "select", width = 500, align = 'start') #dropdown menu to select from

def select_callback(attr, old,new):
    if new: #when new value is selected
        if old: #if old value was not blank - change only plots
            p = p_abs_list[new] #initiate plot
            if new == 'Kopā': #init a total dataframe
                df = combined_data
                df = df.rename(columns = rename_dict)
            else: #init dataframe based on selected value
                df = combined_data.loc[combined_data['testing_lab'] == new]
                df = df.rename(columns = rename_dict)
            data_table = DataTable(columns=Columns, source=ColumnDataSource(df), autosize_mode="fit_columns", aspect_ratio = 1, align = 'start', name = 'table')
            curdoc().get_model_by_name('row').children[0] = p #replace plot
            curdoc().get_model_by_name('row').children[1].children[0] = data_table #replace table
            lg_r_list = {item.renderers[0].name:item.renderers[0] for item in curdoc().get_model_by_name('legend').items}  #filling list of legend renderers after new entry was selected
            for key in lg_r_list.keys():lg_r_list[key].on_change('visible', legend_callback)
        else: #if old value is blank - replace the layout 
            button_1, button_2 = button_dict[new][0], button_dict[new][1] #initiate hide/show buttons
            checkbox = button_dict[new][2]
            button_3 = button_dict[new][3]
            r_slider = button_dict[new][4]
            
            p = p_abs_list[new] #initiate plot
            if new == 'Kopā':
                df = combined_data
                df = df.rename(columns = rename_dict)
            else:
                df = combined_data.loc[combined_data['testing_lab'] == new]
                df = df.rename(columns = rename_dict)
            data_table = DataTable(columns=Columns, source=ColumnDataSource(df), autosize_mode="fit_columns", aspect_ratio = 1, align = 'start', name = 'table')
            c_t_i = column(select,row(button_1, button_2, name = 'buttons'), row(checkbox, r_slider), row(p, Column(data_table, button_3), name = 'row'), name = 'root') #combine controls & plot in a layout
            bokeh_doc.clear() #clear old layout
            bokeh_doc.add_root(c_t_i) #show new layout
            bokeh_doc.get_model_by_name('checkbox').active = []
            lg_r_list = {item.renderers[0].name:item.renderers[0] for item in curdoc().get_model_by_name('legend').items}  #filling list of legend renderers after new entry was selected
            for key in lg_r_list.keys():lg_r_list[key].on_change('visible', legend_callback)
    elif not new:
        c_t_i = column(select, empty_fig) #create a layout with a placeholder column
        bokeh_doc.clear() #clear old layout
        bokeh_doc.add_root(c_t_i) #show new layout
select.on_change('value', select_callback) #map callback to the dropdown object

bokeh_doc.add_root(Column(select, empty_fig, name = 'root')) #add dropdown menu & placeholder figure to the served doc
layouts = bokeh_doc.get_model_by_name('root').children #save layout objects to use later when providing interactivity

button_dict = {} #links buttons to the plots
legend = unq_dict['lineage'] #used to show/hide values on bar plots with buttons

#initating button objects
button_1 = Button(label="Paslēpt visu", button_type="success", width = 250,  align = 'start', name = 'button_1')
button_2 = Button(label="Parādīt visu", button_type="success", width = 250,  align = 'start', name = 'button_2')
checkbox = CheckboxGroup(labels=['Rādīt %'], name = 'checkbox',  align = 'end')
button_3 = Button(label="Eksportēt datus", button_type="success", width = 250,  align = 'center', name = 'button_3')
r_slider = DateRangeSlider(
    value = (date(int(unq_dict['sampling_date'][0][0:4]), int(unq_dict['sampling_date'][0][5:7]), int(unq_dict['sampling_date'][0][-2:])
        ), date(int(unq_dict['sampling_date'][-1][0:4]), int(unq_dict['sampling_date'][-1][5:7]), int(unq_dict['sampling_date'][-1][-2:]))),
    start=date(int(unq_dict['sampling_date'][0][0:4]), int(unq_dict['sampling_date'][0][5:7]), int(unq_dict['sampling_date'][0][-2:])), 
    end=date(int(unq_dict['sampling_date'][-1][0:4]), int(unq_dict['sampling_date'][-1][5:7]), int(unq_dict['sampling_date'][-1][-2:])),
    name = 'slider',
    width = 1000,
    align = 'start')



#defining callback functions
def button_1_callback():
    """Hide all on button press if not already hidden."""
    lg_r_list = {item.renderers[0].name:item.renderers[0] for item in curdoc().get_model_by_name('legend').items}  #filling list of legend renderers after new entry was selected
    df = combined_data.loc[combined_data['lineage']=='Empty']
    df = df.rename(columns = rename_dict)
    for key in lg_r_list.keys():lg_r_list[key].remove_on_change('visible', legend_callback)
    for lg in legend: 
        try: bokeh_doc.get_model_by_name(lg).visible = False
        except: continue
    data_table = DataTable(columns=Columns, source=ColumnDataSource(df), autosize_mode="fit_columns", aspect_ratio = 1, align = 'start', name = 'table')
    curdoc().get_model_by_name('row').children[1].children[0] = data_table #replace table
    for key in lg_r_list.keys():lg_r_list[key].on_change('visible', legend_callback)
    
def button_2_callback():
    """Show all on button press if not already shown."""
    lg_r_list = {item.renderers[0].name:item.renderers[0] for item in curdoc().get_model_by_name('legend').items}  #filling list of legend renderers after new entry was selected
    range_start = datetime.fromtimestamp(r_slider.value[0] / 1e3) #get datetime object from slider current position - range start
    range_end = datetime.fromtimestamp(r_slider.value[1] / 1e3) #get datetime object from slider current position - range end
    list_of_dates = date_range(range_start,range_end)
    
    if curdoc().get_model_by_name('select').value in list(combined_data['testing_lab']):
        df = combined_data.loc[(combined_data['testing_lab'] == curdoc().get_model_by_name('select').value) & (combined_data['sampling_date'].isin(list_of_dates))]
    else: df = combined_data.loc[(combined_data['sampling_date'].isin(list_of_dates))]
    df = df.rename(columns = rename_dict)
    for key in lg_r_list.keys():lg_r_list[key].remove_on_change('visible', legend_callback)
    for lg in legend: 
        try: bokeh_doc.get_model_by_name(lg).visible = True
        except: continue
    data_table = DataTable(columns=Columns, source=ColumnDataSource(df), autosize_mode="fit_columns", aspect_ratio = 1, align = 'start', name = 'table')
    curdoc().get_model_by_name('row').children[1].children[0] = data_table #replace table
    for key in lg_r_list.keys():lg_r_list[key].on_change('visible', legend_callback)

def button_3_callback():
    '''Generate data on the server on-click to be served by other protocol'''
    data_df = pd.DataFrame(curdoc().get_model_by_name('table').source.data)
    data_df.to_csv('/home/user/Desktop/RAKUS/31052021/report.csv',index = False, header=True)
    url = '/home/user/Desktop/RAKUS/31052021/report.csv'
    subprocess.check_call(['libreoffice', url])

    #uncomment when other options are available
    # url = f'{full_path}/4_reports/report_2405/found_samples_20210524-123754.csv'
    # response = requests.get(url)
    # with open('/home/user/Git/cov_flow/downstream/dstr_subprocess/found_samples_20210511-163630.csv', 'r+') as f:
    #     writer = csv.writer(f)
    #     for line in .iter_lines():
    #         writer.writerow(line.decode('utf-8').split(','))

def checkbox_callback(new):
    '''Show relative values for given selection if checked.'''

    range_start = datetime.fromtimestamp(r_slider.value[0] / 1e3) #get datetime object from slider current position - range start
    range_end = datetime.fromtimestamp(r_slider.value[1] / 1e3) #get datetime object from slider current position - range end
    try: start_index = unq_dict['sampling_date'].index(range_start.strftime('%Y-%m-%d')) #if date is covered by data - use the date directly
    except: 
        start_index = unq_dict['sampling_date'].index(nearest(datetime_s_dates,range_start).strftime('%Y-%m-%d')) #if date is not covered by data - switch to nearest date covered
    try: end_index = unq_dict['sampling_date'].index(range_end.strftime('%Y-%m-%d'))+1 #if date is covered by data - use the date directly
    except:
        end_index = unq_dict['sampling_date'].index(nearest(datetime_s_dates,range_end).strftime('%Y-%m-%d')) #if date is not covered by data - switch to nearest date covered
    new_date_range = {'sampling_date':unq_dict['sampling_date'][start_index:end_index]} #new date range to use in plotting
    #generate new plot data - only for selected dates
    if len(bokeh_doc.get_model_by_name('checkbox').active) != 0:
        units = 'rel'
        new_data_dict = {}
        for key in by_lab_lin_dict_rel.keys():
            new_data_dict[key] = {lin:by_lab_lin_dict_rel[key][lin][start_index:end_index] for lin in by_lab_lin_dict_rel[key].keys()} #new data to use in plotting
            for lin in by_lab_lin_dict_rel[key].keys(): #updating legend - removing entries that have no data in selected interval
                if lin != "sampling_date":
                    if sum(by_lab_lin_dict_rel[key][lin][start_index:end_index]) == 0:
                        del new_data_dict[key][lin]
    else:
        units = 'abs'
        new_data_dict = {}
        for key in by_lab_lin_dict.keys():
            new_data_dict[key] = {lin:by_lab_lin_dict[key][lin][start_index:end_index] for lin in by_lab_lin_dict[key].keys()} #new data to use in plotting
            for lin in by_lab_lin_dict_rel[key].keys(): #updating legend - removing entries that have no data in selected interval
                if lin != "sampling_date":
                    if sum(by_lab_lin_dict_rel[key][lin][start_index:end_index]) == 0:
                        del new_data_dict[key][lin]
    #generate new plot for selected dates
    lab = curdoc().get_model_by_name('select').value #get current lab from selection tool
    p = stacked_bar_plot(new_data_dict[lab], new_date_range,'Celms', units=units) #create a stacked bar plot based on rangeslider selection
    #replace the old plot with new one
    curdoc().get_model_by_name('row').children[0] = p #update bokeh document with new plot
    lg_r_list = {item.renderers[0].name:item.renderers[0] for item in curdoc().get_model_by_name('legend').items} #filling list of legend renderers after new interval was selected
    for key in lg_r_list.keys(): lg_r_list[key].on_change('visible', legend_callback)    
    #Get all consecutive dates between two dates chosen on slider
    list_of_dates = date_range(range_start,range_end)
    #Get all rows that contain dates within list_of_dates and match the selected lab
    if curdoc().get_model_by_name('select').value == 'Kopā':
            df = combined_data.loc[combined_data['sampling_date'].isin(list_of_dates)]
            df = df.rename(columns = rename_dict)
    else:
        df = combined_data.loc[(combined_data['testing_lab'] == curdoc().get_model_by_name('select').value) & (combined_data['sampling_date'].isin(list_of_dates))]
        df = df.rename(columns = rename_dict) #rename columns to lv
    data_table = DataTable(columns=Columns, source=ColumnDataSource(df), autosize_mode="fit_columns", aspect_ratio = 1, align = 'end', name = 'table') #generate bokeh table
    curdoc().get_model_by_name('row').children[1].children[0] = data_table #update bokeh table in the layout
    

##custom function to use in rangeslider callback - ref: https://stackoverflow.com/questions/32237862/find-the-closest-date-to-a-given-date
def nearest(items, pivot):
    '''Returns nearest date from list of dates (datetime objects) when passing a datetime object to avoid error when slider reaches dates not covered by data '''
    return min(items, key=lambda x: abs(x - pivot))
datetime_s_dates = [datetime.strptime(date, '%Y-%m-%d') for date in unq_dict['sampling_date']] #to use with nearest()

##custom function to update table values based on rangeslider values - ref: https://www.codegrepper.com/code-examples/python/get+all+dates+between+two+dates+python
def date_range(start, end):
    delta = end - start  # as timedelta
    days = [(start + timedelta(days=i)).strftime('%Y-%m-%d') for i in range(delta.days + 2)] #+2 to include the end date
    return days


def slider_callback(attr, old, new):    
    range_start = datetime.fromtimestamp(r_slider.value[0] / 1e3) #get datetime object from slider current position - range start
    range_end = datetime.fromtimestamp(r_slider.value[1] / 1e3) #get datetime object from slider current position - range end
    try: start_index = unq_dict['sampling_date'].index(range_start.strftime('%Y-%m-%d')) #if date is covered by data - use the date directly
    except: 
        start_index = unq_dict['sampling_date'].index(nearest(datetime_s_dates,range_start).strftime('%Y-%m-%d')) #if date is not covered by data - switch to nearest date covered
    try: end_index = unq_dict['sampling_date'].index(range_end.strftime('%Y-%m-%d'))+1 #if date is covered by data - use the date directly
    except:
        end_index = unq_dict['sampling_date'].index(nearest(datetime_s_dates,range_end).strftime('%Y-%m-%d')) #if date is not covered by data - switch to nearest date covered
    new_date_range = {'sampling_date':unq_dict['sampling_date'][start_index:end_index]} #new date range to use in plotting
    #generate new plot data - only for selected dates
    if len(bokeh_doc.get_model_by_name('checkbox').active) != 0:
        units = 'rel'
        new_data_dict = {}
        for key in by_lab_lin_dict_rel.keys(): #updating legend - removing entries that have no data in selected interval
            new_data_dict[key] = {lin:by_lab_lin_dict_rel[key][lin][start_index:end_index] for lin in by_lab_lin_dict_rel[key].keys()} #new data to use in plotting
            for lin in by_lab_lin_dict_rel[key].keys(): 
                if lin != "sampling_date":
                    if sum(by_lab_lin_dict_rel[key][lin][start_index:end_index]) == 0:
                        del new_data_dict[key][lin]
    else:
        units = 'abs'
        new_data_dict = {}
        for key in by_lab_lin_dict.keys(): #updating legend - removing entries that have no data in selected interval
            new_data_dict[key] = {lin:by_lab_lin_dict[key][lin][start_index:end_index] for lin in by_lab_lin_dict[key].keys()} #new data to use in plotting
            for lin in by_lab_lin_dict_rel[key].keys(): 
                if lin != "sampling_date":
                    if sum(by_lab_lin_dict_rel[key][lin][start_index:end_index]) == 0:
                        del new_data_dict[key][lin]
    #generate new plot for selected dates
    lab = curdoc().get_model_by_name('select').value #get current lab from selection tool
    p = stacked_bar_plot(new_data_dict[lab], new_date_range,'Celms', units=units) #create a stacked bar plot based on rangeslider selection
    #replace the old plot with new one
    curdoc().get_model_by_name('row').children[0] = p #update bokeh document with new plot
    lg_r_list = {item.renderers[0].name:item.renderers[0] for item in curdoc().get_model_by_name('legend').items} #filling list of legend renderers after new interval was selected
    for key in lg_r_list.keys(): lg_r_list[key].on_change('visible', legend_callback)                        
    #Get all consecutive dates between two dates chosen on slider
    list_of_dates = date_range(range_start,range_end)
    #Get all rows that contain dates within list_of_dates and match the selected lab
    if curdoc().get_model_by_name('select').value == 'Kopā':
            df = combined_data.loc[combined_data['sampling_date'].isin(list_of_dates)]
            df = df.rename(columns = rename_dict)
    else:
        df = combined_data.loc[(combined_data['testing_lab'] == curdoc().get_model_by_name('select').value) & (combined_data['sampling_date'].isin(list_of_dates))]
        df = df.rename(columns = rename_dict) #rename columns to lv
    data_table = DataTable(columns=Columns, source=ColumnDataSource(df), autosize_mode="fit_columns", aspect_ratio = 1, align = 'end', name = 'table') #generate bokeh table
    curdoc().get_model_by_name('row').children[1].children[0] = data_table #update bokeh table in the layout

#mapping callbacks to button objects
button_1.on_click(button_1_callback)
button_2.on_click(button_2_callback)
button_3.on_click(button_3_callback)
checkbox.on_click(checkbox_callback)
r_slider.on_change("value_throttled", slider_callback) #throttled allows to avoid updating upon each slider movement

#mapping buttons to data
for key in by_lab_lin_dict_rel.keys(): button_dict[key] = [button_1, button_2, checkbox, button_3, r_slider]
