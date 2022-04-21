#!/bin/bash
#PBS -n read_library
#PBS -l walltime=00:05:30
#PBS -l procs=3
#PBS -q batch
#PBS -j oe
#PBS -A rakus

# THE SCRIPT IS USED TO GENERATE LIBRARY COMPOSITION REPORTS FROM DEMULTIPLEXED, UNTRIMMED FASTQ FILE.

module load singularity
fastq_screen_sif_path=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "fastq_screen.sif")
file_path=$(sed 's/ /\\ /g' <<< ${1}) #PATH TO READ_1 FASTQ FILE (ESCAPING SPACES IN FILE NAMES)
base_dir=${2} #PATH TO WORKING DIRECTORY WHERE TEMPORARY AND FINAL OUTPUT FILES ARE SAVED
now=$(date +"%m_%d_%Y") #TO BE USED LATER IN NAMED REPORT FOLDERS
output_path=$(find ${base_dir} -maxdepth 2 -type d -name "fastq_screen_output")
output_path=${output_path#"/HOME/GROUPS/NMRL/"} #CONVERT ABSOLUTE PATH TO RELATIVE
cd /home/groups/nmrl/

#CALCULATE DOWNSAMPLING SIZE - 20% OF THE ORIGINAL DATASET (4*5 IN THE DOWNSAMPLE FORMULA DENOMINATOR) - 4 - USED TO COUNT READS IN FASTQ.GZ FILE (NEW READ ON EACH 4TH ROW); 5 - TAKING 20% OF READ COUNT
downsample=$(zcat $file_path | echo $((`wc -l`/20)))
singularity run $fastq_screen_sif_path fastq_screen --subset $downsample -conf db/db-fastq-screen/fastq_screen.conf $file_path --outdir $output_path/
