#!/bin/bash

report_path=${1}
data_path=${2}
regex=$(cat ./current_regex.txt)
ls ${report_path} | grep -oP ${regex} > ${report_path}/processed_samples.txt
awk '{print $0".*.fastq.gz"}' ${report_path}/processed_samples.txt > ${report_path}/filter_list.txt
cd ${data_path}
ls -R | grep -f ${report_path}filter_list.txt > ${report_path}fastq_list.txt
rm ${report_path}processed_samples.txt
rm ${report_path}filter_list.txt
