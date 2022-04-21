#!/bin/sh
# This invokes screen if not already inside a screen which invokes /bin/bash, which invokes the script again.
fastq_folder_path=${1}
output_dir=${2}
#if [ -z "$STY" ]; then exec screen -dm -S c19_processing /bin/bash "$0"; fi
python3 /home/groups/nmrl/cov_analysis/fastq_processing/run_fastq_screen.py -f ${fastq_folder_path} -n 8 -o ${output_dir} > ${output_dir}/fastq_screen.log &