import pandas as pd, sys, argparse

#THE SCRIPT IS USED TO GENERATE BATCH SUMMARY FILE TO BE MERGED WITH DATABASE FILE

parser = argparse.ArgumentParser(description='A script to generate batch report from new mutation report csv file.')
parser.add_argument('-p', '--pipe', metavar='\b', help = 'Full path to the mutation_report.csv', default=None, required=False)

#IF NO ARGS PROVIDED - PRINT HELP & EXIT
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

#PARSING ARGUMENTS
pipe_report = args.pipe #PATH TO PIPELINE REPORT

#SETTING OUTPUT PATH BASED ON INPUT (PIPE_REPORT)
output_path = '/'.join(pipe_report.split('/')[:-1])

#READ DATA
pipe_data = pd.read_csv(pipe_report)

#DEFINE DATABASE COLUMNS (REPLACE BY READING COLUMNS FROM SUMMARY?)
columns = [
    'receiving_lab_sample_id',
    'processing_id',
    'seq_institution',
    'testing_lab_sample_id',
    'lineage',
    'genome_length',
    'genome_N_percentage',
    'genome_GC_content',
    'sampling_date',
    'testing_lab',
    'age',
    'gender',
    'district',
    'sample_type', 
    'AVERAGE_COVERAGE',
    'MEDIAN_COVERAGE',
    'MAPPED_FRACTION',
    'READS_MAPPED',
    'TOTAL_READS',
    'seq_date'
]

#EXTRACT METRICS FROM MUTATION REPORT
pipe_data['testing_lab_sample_id'] = pipe_data['receiving_lab_sample_id'] #TO MATCH COLUMN USED TO STORE DOUBLE-LABELLED IDS (ARTIFACT)
pipe_data['sample_type'] = [0 for _ in range(len(pipe_data))] #TO MATCH COLUMN STORING SAMPLE TYPE (ARTIFACT)
batch_report = pipe_data[columns]

#FINALIZING THE DATAFRAME
batch_report.fillna(0, inplace=True) #STANDARDIZE EMPTY RECORDS
batch_report.drop_duplicates(subset="receiving_lab_sample_id", keep='first', inplace=True) #REMOVE RECORDS WHERE SAMPLE ID IS DUPLICATED
batch_report.to_csv(f'{output_path}/batch_report.csv', header=True, index=False) #SAVE TO CSV
