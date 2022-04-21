###SCOPE
"""The script is used to run SARS-CoV2 pipeline on HPC."""

###IMPORTS
import os, subprocess, re, time, sys, argparse, pandas as pd
from collections import deque
from pathlib import Path


###ARGUMENT PARSING CONFIGURATION
parser = argparse.ArgumentParser(description='A script to run the analysis of SARS-CoV19 sequencing data in fastq format.') #ARGPARSER OBJECT TO PROVIDE COMMAND-LINE FUNCTIONALITY
req_arg_grp = parser.add_argument_group('required arguments') #OPTIONS WILL BE DISPLAYED AS SEPARATED REQUIRED ARG GROUP IN HELP MESSAGE


###ADDING COMMAND-LINE ARGUMENTS
req_arg_grp.add_argument('-f', '--fastq', metavar='\b', help = 'Full path to the folder containing fastq files to process. All files from folder and subfolders will be included in the analysis.', required=True)
###TO BE REFACTORED: parser.add_argument('-p', '--skip_processed', metavar='\b', help = 'Path to list of processed fastq files (use find_processed.sh to generate list).', default=False, required=False)
parser.add_argument('-n', '--n_samples', metavar='\b', help = 'Number of samples to be processed in-parallel on different nodes of the cluster (8 by default)', default=8, required=False)
req_arg_grp.add_argument('-s', '--summary', metavar='\b', help = 'Path to summary file containing latest processing id.', default=None, required=True)


###IF NO ARGUMENTS PROVIDED - SHOW HELP & EXIT
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)


###PARSE ARGUMENTS
args = parser.parse_args()


###PATHS & CONFIGS
full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[0] #PATH TO SCRIPT LOCATION
tool_path = f'{full_path}/fastq_processing/tools/' #SOFTWARE STORED OUTSIDE CONTAINER (SNPEFF & ABRA JARS)
fastq_path = args.fastq #PATH TO FOLDER WITH FASTQ FILES TO PROCESS (PARENT DIR)
report_time = time.strftime("%Y_%m_%d") #TIMESTAMP
sequencing_lab = f"NMRL_ECDC" #DEFALT SEQ LAB TO USE FOR CONSISTANT FOLDER NAMING
sequencing_time = str(None) #DEFAULT SEQ TIME
home_path = os.path.expanduser('~/') #PATH TO USER HOME DIR


###ADJUSTING SEQ_LAB FOR NMRL SAMPLES
if "NMRL" in fastq_path: 
    sequencing_lab = f"NMRL_LIC"


###ADJUSTING SEQ_DATE FOR NMRL SAMPLES
found_date = re.search(r"[0-9]{4}_[0-9]{2}_[0-9]{2}",fastq_path)
if found_date is not None:
    sequencing_time = found_date.group(0)


###MOVE DATA MOUNT POINT TO USER HOME DIR TO EASE THE PROCESSING AND REDIRECT FASTQ_FILE PATH
os.system(f'mv -n {fastq_path} {home_path}/')
if fastq_path[-1] == "/":
    fastq_path = fastq_path[:-1]
    dir_name = fastq_path.rsplit('/',1)[-1]
    fastq_path = f'{home_path}{dir_name}'
else:
    dir_name = fastq_path.rsplit('/',1)[-1]
    fastq_path = f'{home_path}{dir_name}'

###PREPARING SAMPLE IDS AND FILE PATHS TO SERVE TO PIPELINE
read_dict = {}
file_list = []
for (root,dirs,files) in os.walk(fastq_path, topdown=True): #GET LIST OF FASTQ FILE PATHS (FROM PARENT DIR & SUBDIR)
    for name in files: 
        if 'fastq.gz' in name: file_list.append(os.path.join(root, name)) 


###CURRENTLY NOT WORKING
# if args.skip_processed:
#     with open(args.skip_processed, "r+") as processed_file: 
#         total_length = len(file_list)
#         name_list = processed_file.read()
#         for file in file_list:
#             if file.split("/")[-1] in name_list:
#                 print(f'{file} is already processed and is excluded from the analysis')
#                 file_list.remove(file)
#         print(f'{total_length - len(file_list)} processed files were excluded from the analysis.')


###GET LATEST INHOUSE ID FROM SUMMARY FILE
summary_path = args.summary
latest_id = list(pd.read_csv(summary_path)['processing_id'])[-1]


###DIFFERENCE BETWEEN LOCAL AND EUROFIN SAMPLES
if file_list[0].split("/")[-1][0] == "C":
    seq_id_list = set(file.split("/")[-1].split("_",1)[-1].split("_",1)[-1].rsplit("_",2)[0] for file in file_list)
else:
    seq_id_list = set(file.split("/")[-1].split("_",1)[0] for file in file_list)


###GENERATE LIST OF PROCESSING_IDS
processing_id_list = [latest_id[:3]+str(int(latest_id[3:])+i+1).zfill(6) for i in range(len(seq_id_list))]

###GENERATE NEW LIST OF COMBINED ID
seq_id_list_new = {ids[0]:f'{ids[1]}_{ids[0]}' for ids in zip(seq_id_list,processing_id_list)}


###ADD READ_1 AND READ_2 PATH TO LIST MATCHED TO SAMPLE ID
for id in seq_id_list: read_dict[id]=[] #INIT DICT TO STORE PATHS TO READ1 AND READ2 FILES
for id in read_dict.keys(): #ITERATE OVER SAMPLE IDS
    for file in file_list: #ITERATE OVER FILE PATHS IN FILE LIST
        if f'{id}_' in file: #IF FILE PATH CONTAINS ID
            read_dict[id].append(file) #ADD FULL PATH TO FASTQ FILE TO DICT
            if len(read_dict[id]) == 2: #IF PATHS TO BOTH READ 1 & READ 2 FILES ARE IN THE DICT
                break #STOP LOOKING FOR PATHS FOR THIS ID


###REPLACE ORIGINAL ID WITH PROCESSING ID
new_read_dict = {}
for key in read_dict.keys():
    new_read_dict[seq_id_list_new[key]] = read_dict[key]
read_dict = new_read_dict


###FILTERING READS THAT ARE IMPROPERLY PAIRED
key_list = list(read_dict.keys())
for key in key_list:
    if len(read_dict[key]) < 2:
        del read_dict[key]


###SETTING WD TO SUBPROCESS TO USE SUBPROCESSING SCRIPTS
os.chdir(f'{full_path}/fastq_processing/subprocess')


###PREPARING WORKING DIRECTORIES FOR NEW USER
if not os.path.isdir(f'{home_path}/output_final/'):
    os.mkdir(f'{home_path}/output_final/')
if not os.path.isdir(f'{home_path}/output_temp/'):
    os.mkdir(f'{home_path}/output_temp/')
if not os.path.isdir(f'{home_path}/output_final/{sequencing_lab}-{sequencing_time}-{report_time}'):
    os.mkdir(f'{home_path}/output_final/{sequencing_lab}-{sequencing_time}-{report_time}')


###START FASTQSCREEN
subprocess.check_call(['../screen_run_fastqscreen.sh', f'{fastq_path.split("/")[-1]+"/"}', f'{home_path}/output_final/{sequencing_lab}-{sequencing_time}-{report_time}'])

###CURRENTLY NOT WORKING
# if args.skip_processed:
#     updated_dict = {}
#     for id in read_dict.keys():
#         if len(read_dict[id]) ==2: updated_dict[id] = read_dict[id]
#     read_dict = updated_dict


def submit_job(sample_id,que):
    '''Function is used to submit job to hpc in a loop'''
    read_dict[sample_id] = sorted(read_dict[sample_id]) #TO ENSURE THAT READ_1 IS SUPPLIED TO THE ANALYSIS SCRIPT BEFORE READ_2
    read_1 = read_dict[sample_id][0]
    read_2 = read_dict[sample_id][1]
    print(sample_id, read_1, read_2)
    print(sample_id, f'{len(que)}/{len(list(read_dict.keys()))}')
    print(f'{read_1} {read_2} {sequencing_lab}-{sequencing_time}-{report_time} {tool_path} {sample_id} {home_path}')
    subprocess.check_call(['qsub', '-F', f'{read_1} {read_2} {sequencing_lab}-{sequencing_time}-{report_time} {tool_path} {sample_id} {home_path}', 'process_fastq.sh'])


###PYTHON QUE PARAMETERS
sample_limit = args.n_samples #NUMBER OF SAMPLES TO BE PROCESSED AT THE SAME TIME ON DIFFERENT NODES
sleep_time = 30 #SLEEP TIME IN SECONDS TO WAIT BEFORE ADDING NEW SAMPLES TO QUE
sample_que = deque(list(read_dict.keys())) #CREATE A QUE OF SAMPLES


###COPY SCRIPTS TO USER HOME DIR TO BE ABLE TO USE THEM FROM PIPELINE
os.system(f'cp -r {full_path}/fastq_processing/resources {home_path}/')
os.system(f'cp {full_path}/fastq_processing/subprocess/depth_plot.py {home_path}/resources/')
os.system(f'cp {full_path}/fastq_processing/subprocess/vcf_to_csv_cmd.py {home_path}/resources/')


###SAMPLE SUBMISSION LOOP
while len(sample_que) > 0: #UNTIL THERE ARE SAMPLES TO BE PROCESSED
    sample_id = sample_que.pop() #POP RIGHT-MOST SAMPLE ID FROM QUE (QUE IS SHORTENED BY ON ELEMENT)
    while len([temp_dir for temp_dir in os.listdir(f'{home_path}/output_temp/') if 'temp_files' in temp_dir]) >= int(sample_limit): #IF THERE ARE ALREADY {SAMPLE_LIMIT} JOBS IN PROGRESS
        print(f'Job submission halted for {sleep_time} seconds - {sample_limit} samples being processed at the moment.')
        time.sleep(sleep_time) #WAIT FOR {SLEEP_TIME SECONDS}
    os.mkdir(f'{home_path}/output_temp/temp_files_{sample_id}')
    submit_job(sample_id,sample_que)


####IF THE FLOW OF EXCECUTION IS HERE THEN ALL SAMPLES IN THE QUE ARE SUBMITTED FOR PROCESSING (LOOP IS COUNTING FOLDERS IN TEMP DIR)
while len([temp_dir for temp_dir in os.listdir(f'{home_path}/output_temp/') if 'temp_files' in temp_dir]) > 0: #WAITING FOR THE PROCESSING OF FINAL SUBMITTED SAMPLES TO COMPLETE BEFORE STARTING CLEANUP
    print(f'Waiting for the processing to finish - next check in {sleep_time} seconds')
    time.sleep(sleep_time)

###CHECK IF FASTQC SUBPROCESS IS COMPLETE
#create lock file in fastqscreen python script
#when fastqscreen complete - remove lock file (with fastqscreen python script)
#in get_covid - check if lock file still present - if there - wait, else - continue

fscreen_lock = os.path.isfile(f"{home_path}/fastq_screen_lock.file")
while fscreen_lock:
    print(f'Waiting for the fastqscreen to finish - next check in {sleep_time} seconds')
    time.sleep(sleep_time)
    fscreen_lock = os.path.isfile(f"{home_path}/fastq_screen_lock.file")

###RUNNING CLEANUP TASKS
os.system(f"mv -n {full_path}/fastq_processing/subprocess/find_covid19_mutations* {home_path}/output_final/{sequencing_lab}-{sequencing_time}-{report_time}/") #move log files to the output folder
os.system(f'rm -fr {home_path}/resources')
os.system(f'chmod -R 775 ~/output_final/{sequencing_lab}-{sequencing_time}-{report_time}/')
os.system(f'mv -n {home_path}/output_final/{sequencing_lab}-{sequencing_time}-{report_time} {full_path}/covid_output/')
os.system(f'mv -n {fastq_path} {full_path}/covid_input/')
