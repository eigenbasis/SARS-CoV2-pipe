import sys, os, pandas as pd, re, time, datetime, concurrent.futures, subprocess
from datetime import datetime
from pathlib import Path
import comb_fasta_stat_calc as get_stats

'''Paths&setups'''
timestr_fasta = time.strftime("%m_%d_%Y")
full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[0] 
f_value = 65 #percentage of reads that support a mutation at given position (%)
folder_path = f'{full_path}/reports/report_{timestr_fasta}' 
subprocess_path = f'{full_path}/subprocess'
metadata_path = f'{full_path}/resources/metadata'


'''Controllers'''
skip_excel_mining = False 
run_pangolin = True
update_pangolin = False
skip_sequence_stats = False
rename_fasta_header = True
metadata_file_path= sys.argv[1]
if metadata_file_path == 'None':  skip_excel_mining = True

'''Metadata'''
if not skip_excel_mining: #if metadata path is provided, extract metadata
    all_excel_dump_df = pd.read_csv(metadata_file_path)
else: all_excel_dump_df = None

'''Running pangolin'''
if run_pangolin: #providing correct path for run_pangolin.sh & using subprocess to run it with options from control variables
    if update_pangolin: os.chdir(subprocess_path), subprocess.check_call(["./run_pangolin.sh", str(1), str(1)]) 
    elif not update_pangolin: os.chdir(subprocess_path), subprocess.check_call(["./run_pangolin.sh", str(1), str(0)]) 

elif not run_pangolin: #if pangolin typing option skipped - check if pangolin udpate is required
    if update_pangolin: os.chdir(subprocess_path), subprocess.check_call(["./run_pangolin.sh", str(0), str(1)]) 
    elif not update_pangolin: os.chdir(subprocess_path), subprocess.check_call(["./run_pangolin.sh", str(0), str(0)])

if run_pangolin: #if running pangolin - extract path to combined fasta file & pango report path
    combined_path = [path for path in os.listdir(f'{full_path}/reports/report_{timestr_fasta}') if 'combined.fasta' in path][0]
    combined_fasta_path = f'{full_path}/reports/report_{timestr_fasta}/{combined_path}'
    pango_report_path = f'{full_path}/reports/report_{timestr_fasta}/{timestr_fasta}_lineage_report.csv'

'''Sequence stats'''
if not skip_sequence_stats: sequence_stats_df = get_stats.get_nt_counting_stats(combined_fasta_path, skip_sequence_stats, rename_fasta_header) #generating sequence statistics with option to skip and process >LVA000_2104047246 and >hCoV-19/Latvia/2103098474/2021 type fasta headers

'''Filters'''
file_paths = [f'{folder_path}/{file_name}' for file_name in os.listdir(folder_path) if 'ann.csv' in file_name] #list of mutation csv reports for each processed sample
coverage_paths = [f'{folder_path}/{file_name}' for file_name in os.listdir(folder_path) if 'seq_depth.txt' in file_name] #list of coverage depth reports for each processed sample
mapped_paths = [f'{folder_path}/{file_name}' for file_name in os.listdir(folder_path) if 'mapped_report.txt' in file_name] #list of samtools flagstats report for each processed sample

with open(f"{Path(os.path.dirname(os.path.realpath(__file__))).parents[0]}/resources/report_filters.txt", "r+") as filter_file:
    filter_list = {filter.split(" ")[0]:set(filter.split(" ")[1].split(",")) for filter in filter_file.read().split("\n")[:-1]} #saves filter_name:filter_set pair to a dictionary

'''Functions'''
def csv_info_extractor(file_path, filter_list, f_value):
    """The function is used to exctract mutation data from vcf-based csv reports based on set of provided filters"""

    sample_id = re.search(r'([0-9]{10}|[0-9]{9})',file_path).group(0)    
    df = pd.read_csv(file_path)
    result_row = {"SAMPLE_ID":sample_id}

    for key in filter_list.keys(): 
        mut_df = df.loc[(df["AMINO_ACID_CHANGE"].isin(filter_list[key])) & (df["FREQUENCY"] >= f_value)]
        if len(mut_df["P_ERR_MUT_CALL"]) == len(filter_list[key]): #if the number of extracted (unique) rows matches the number of mutation_names for a given filter_name
            result_row[key] = str(1) #assume that the mutation specified by the filter was found - add filter_name:1 pair to the result_row dict
        else: result_row[key] = str(round(len(mut_df["P_ERR_MUT_CALL"])/len(filter_list[key]), 2)) #if the number of extracted (unique) rows does not match the number of mutation_names for a given filter_name - add filter_name:%match pair to the result_row dict
    
    return result_row

def depth_info_extractor(file_path):
    """The function is used to extract average coverage depth from samtools depth report"""

    data = pd.read_csv(file_path, delimiter='\t').iloc[:,2] #convert depth report to a dataframe
    sample_id = re.search(r'([0-9]{10}|[0-9]{9})',file_path).group(0) #extract sample if from file name
    result_row = {'SAMPLE_ID':sample_id, 'AVERAGE_COVERAGE':round(sum(data)/len(data),2)} #map average coverage to sample id in a dict
    
    return result_row

def mapped_info_extractor(file_path):
    """The function is used to extract mapped read data from samtools flagstats report file."""

    with open(file_path, "r+") as file: data = file.readlines() #read flatstats file
    sample_id = re.search(r'([0-9]{10}|[0-9]{9})',file_path).group(0) #extract sample id from file name
    result_row = {
        'SAMPLE_ID':sample_id, 
        'TOTAL_READS':data[0].split(" ")[0], 
        'READS_MAPPED':data[4].split(" ")[0], 
        'MAPPED_FRACTION':round(int(data[4].split(" ")[0])/int(data[0].split(" ")[0]),2)} #map samtools flagstats output to sample id in a dict
    
    return result_row

def mutation_report_generator(file_paths, coverage_paths, output_path, f_value, all_excel_dump_df=None, sequence_stats_df=None, pango_report_path=None):
    """The function is used to generate mutation report based on data extracted from mutations reports for individual samples."""

    result_df = pd.DataFrame() #init empty dataframe for the report
    
    with concurrent.futures.ProcessPoolExecutor() as executor: #applying the extractor function in-parallel on different cores
        start_time = time.time() #to estimate the time it took to process the files
        results = [executor.submit(csv_info_extractor, file_path, filter_list, f_value) for file_path in file_paths] #submitting function calls to different processes

        for f in concurrent.futures.as_completed(results): #collecting processing results to the dataframe
            print(f'{f.result()["SAMPLE_ID"]} - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.now())}') #printing the time when the processing was finished for a given sample
            result_df = result_df.append(f.result(), ignore_index=True) #add the result_row to the mutation report
        print(f"Processed {len(file_paths)} mutation reports in --- %s seconds ---" % round(time.time() - start_time, 2)) #the time it took to process all the mutation report files from file_path_set    
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        start_time = time.time() #to estimate the time it took to process the files
        results = [executor.submit(depth_info_extractor, file_path) for file_path in coverage_paths] #submitting function calls to different processes
        coverage_df = pd.DataFrame() #init empty dataframe to store samtools flagstats data for multiple samples

        for f in concurrent.futures.as_completed(results): #collecting processing results to the dataframe
            print(f'{f.result()["SAMPLE_ID"]} - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.now())}') #printing the time when the processing was finished for a given sample
            coverage_df = coverage_df.append(f.result(), ignore_index=True) #add the result_row to the report df
        print(f"Processed {len(file_paths)} coverage reports in --- %s seconds ---" % round(time.time() - start_time, 2))

    result_df = pd.merge(result_df.applymap(str), coverage_df.applymap(str), how="left", on="SAMPLE_ID") #add coverage depth data to the mutation report

    with concurrent.futures.ProcessPoolExecutor() as executor:
        start_time = time.time() #to estimate the time it took to process the files
        results = [executor.submit(mapped_info_extractor, file_path) for file_path in mapped_paths] #submitting function calls to different processes
        mapped_df = pd.DataFrame()

        for f in concurrent.futures.as_completed(results): #collecting processing results to the dataframe
            print(f'{f.result()["SAMPLE_ID"]} - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.now())}') #printing the time when the processing was finished for a given sample
            mapped_df = mapped_df.append(f.result(), ignore_index=True) #add the result_row to the report df
        print(f"Processed {len(file_paths)} mapping reports in --- %s seconds ---" % round(time.time() - start_time, 2))

    result_df = pd.merge(result_df.applymap(str), mapped_df.applymap(str), how="left", on="SAMPLE_ID") #add samtools flagstats data to the mutation report

    column_reorder_list = ['SAMPLE_ID'] + list(filter_list.keys()) + ['AVERAGE_COVERAGE', 'MAPPED_FRACTION', 'READS_MAPPED', 'TOTAL_READS'] #specifying the column order in the mutation report
    result_df = result_df[column_reorder_list] #reordering columns based in order of filters in report_filters.txt file
    result_df.rename(columns={'SAMPLE_ID':"receiving_lab_sample_id"}, inplace=True) #for merging purposes

    if all_excel_dump_df is not None: #if metadata mining option selected - add metadata to the mutation report
        result_df = pd.merge(result_df.applymap(str), all_excel_dump_df.applymap(str), how="left", on="receiving_lab_sample_id")

    if sequence_stats_df is not None: #if statistics calculation option selected - add sequence statistics data to the mutation report
        result_df = pd.merge(result_df.applymap(str), sequence_stats_df.applymap(str), how="left", on="receiving_lab_sample_id")

    if pango_report_path is not None: #if pangolin typing option selected - add pangolin lineage data to the mutation report
        pango_report_df= pd.read_csv(pango_report_path) #read pangolin report to a dataframe

        if re.match(r'[A-Z]{2}[0-9]{4}\.B', pango_report_df['taxon'][0]): pango_id_dict = {taxon:taxon for taxon in pango_report_df['taxon']}#map taxot to sample id in a dict
        elif rename_fasta_header: pango_id_dict = {taxon:taxon.split("/")[-2] for taxon in pango_report_df["taxon"]} #option for GISAID-accepted headers
        else: pango_id_dict = {taxon:taxon.split("_")[1] for taxon in pango_report_df["taxon"]} #map taxot to sample id in a dict if fasta headers were not changed

        pango_report_df["receiving_lab_sample_id"] = pango_report_df["taxon"].map(pango_id_dict) #adding ([0-9]{9}|[0-9]{10}) format sample ids to the pango report dataframe
        pango_report_df = pango_report_df[["receiving_lab_sample_id", 'lineage']] #keep only lineage matched against sample ids
        result_df = pd.merge(result_df.applymap(str), pango_report_df.applymap(str), how="left", on="receiving_lab_sample_id") #adding pangolin report data to the mutation report
    
    result_df.to_csv(f'{output_path}/{"{0:%Y-%m-%d-%H-%M}".format(datetime.now())}_ft_{f_value}_mutation_report.csv', header=True, index=False) #generating the report

mutation_report_generator(file_paths, coverage_paths, folder_path, f_value, all_excel_dump_df, sequence_stats_df, pango_report_path)