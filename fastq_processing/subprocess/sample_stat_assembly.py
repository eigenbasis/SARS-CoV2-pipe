import pandas as pd, os, sys
from pathlib import Path

'''The script is used to automatically update summary table that is used for data visualization each time new batch of samples is processed with a metadata file.'''
full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[1] #one-path up path
summary_path = f'/mnt/c/cov_seq/4_reports/summary/sample_summary_stats1.csv' #static path to the base statistics file
# summary_path = f'{full_path}/test_summary_stats.csv' #path to summary file
report_path = sys.argv[1] #path to freshly generated report containing metadata is expected
summary_df, report_df = pd.read_csv(summary_path), pd.read_csv(report_path) #read both summary table and new report

for col in report_df.columns: #add new column into summary table if it is not there but found in new report and full it with 0-s
    if col not in summary_df.columns:
        summary_df[col] = ['0' for i in range(len(summary_df))]

summary_df = summary_df.append(report_df) #append new report data to the summary table
summary_df.drop_duplicates(inplace=True) #remove duplicated values from the summary table
summary_df[sorted(summary_df.columns)] #sort columns in summary table alphabetically
summary_df.to_csv(summary_path, header=True, index=False) #replace olde summary table with the updated summary table
