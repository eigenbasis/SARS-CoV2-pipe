import os, subprocess, re, time, sys, argparse, pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(description='A script to run the analysis of SARS-CoV19 sequencing data in fastq format.') #argparser object to provide command-line functionality
req_arg_grp = parser.add_argument_group('required arguments') #add optional argument group
req_arg_grp.add_argument('-f', '--fastq', metavar='\b', help = 'Full path to the folder containing fastq files to process.', required=True) #required fastq file argument
parser.add_argument('-m', '--metadata', metavar='\b', help = 'Name of the metadata file under resources/metadata directory.', default=None, required=False) #optional metadata argument
parser.add_argument('-v', '--visualize', metavar='\b', help = 'True/False if data visualization is required/not required.', default=False, required=False) #optional visualization argument

if len(sys.argv)==1: #if no command-line arguments provided - display help and stop script excecution
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args() #args list from command-line input

full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[0] #fastq processing folder
plot_data = f'{Path(os.path.dirname(os.path.realpath(__file__))).parents[0]}/test_summary_stats.csv' #path to summary file to be used in plotting
# tool_path = f'/mnt/c/cov_seq/tools' #path to tools
tool_path = f'/home/user/tools'
fastq_path = args.fastq #path to fastq files from command-line argument

report_time = time.strftime("%m_%d_%Y")
if args.metadata is not None: metadata_path = f'{full_path}/fastq_processing/resources/metadata/{args.metadata}' #if metadata provided - read metadata file
else: metadata_path = None

read_dict = {}
file_list = []
for (root,dirs,files) in os.walk(fastq_path, topdown=True): # iterate over all subdir under fastq_path and add path to fastq files to list
    for name in files: 
        if 'fastq.gz' in name: file_list.append(os.path.join(root, name)) 

seq_id_list = set(re.search(r'([0-9]{10}|[0-9]{9})', file).group(0) for file in file_list) # get set of unique sample ids form fastq file list
date_dict = {id:'' for id in seq_id_list} #init dict to store sequencing date for each sample
for id in seq_id_list: read_dict[id]=[] #init dict to store paths to read1 and read2 files
for file in file_list: #add read_1 and read_2 path to list matched to sample id
    read_dict[re.search(r'([0-9]{10}|[0-9]{9})', file).group(0)].append(file) 
    if  re.search(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', file):
        date_dict[re.search(r'([0-9]{10}|[0-9]{9})', file).group(0)] = re.search(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', file).group(0) 
    else: 
        date_dict[re.search(r'([0-9]{10}|[0-9]{9})', file).group(0)] = '0'

os.chdir(f'{full_path}/fastq_processing/subprocess')

if not os.path.isdir(f'{full_path}/fastq_processing/reports'): #add reports folders if they do not exist
    os.mkdir(f'{full_path}/fastq_processing/reports')
if not os.path.isdir(f'{full_path}/fastq_processing/reports/report_{report_time}'):
    os.mkdir(f'{full_path}/fastq_processing/reports/report_{report_time}')

for id in list(read_dict.keys()): #running the analysis script
    read_dict[id] = sorted(read_dict[id]) #to ensure that read_1 is supplied to the analysis script before read_2
    read_1 = read_dict[id][0]
    read_2 = read_dict[id][1]
    print(id, read_1, read_2)
    subprocess.check_call(["./process_fastq.sh", read_1, read_2, f'report_{report_time}', tool_path])


os.system(f'python mutation_report.py {metadata_path}') #generate combined mutation report when all samples in sample list were processed
report_path = None
for file in os.listdir(f'{full_path}/fastq_processing/reports/report_{report_time}'): #adding sequencing date data to the combined mutation report
    if 'mutation_report.csv' in file:
        report_path = f'{full_path}/fastq_processing/reports/report_{report_time}/{file}'
        mut_df = pd.read_csv(report_path).astype(str)
        mut_df['seq_date'] = mut_df['receiving_lab_sample_id'].map(date_dict)
        seq_date = mut_df.pop('seq_date')
        mut_df.insert(1, 'seq_date', seq_date)
        mut_df.to_csv(report_path, header=True, index = False)
if report_path is None:
    print('No mutation report found.')
    sys.exit()


if metadata_path is not None: os.system(f'python sample_stat_assembly.py {report_path}') #update summary statistics table from new mutation report

if args.visualize: #run data visualization
    if len(pd.read_csv(plot_data)) >= 2: #if less than 2 records are in the summary stats table - do not show dashboard
        os.chdir(f'{full_path}/fastq_processing/plotting/')
        subprocess.check_call(['bokeh', 'serve', '--show', 'lablin_data_serve.py', '--args', plot_data])
    else:
        print('Data does not contain enough samples to create the dashboard (2 or more are required).')
