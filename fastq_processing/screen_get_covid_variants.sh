#!/bin/sh
# This invokes screen if not already inside a screen which invokes /bin/bash, which invokes the script again.
fastq_folder_path=${1}
#processed_path=${2}
#if [ -z "$STY" ]; then exec screen -dm -S c19_processing /bin/bash "$0"; fi
python3 /home/groups/nmrl/cov_analysis/fastq_processing/get_covid_variants.py -f ${fastq_folder_path} -n 16 -s /home/groups/nmrl/cov_analysis/downstream/summary_file_19_04_2022.csv > ~/output_final/output.log &
