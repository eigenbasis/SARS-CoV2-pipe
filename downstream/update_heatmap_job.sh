#!/bin/bash
#PBS -N update_mut_matrix
#PBS -l walltime=00:30:00
#PBS -l procs=20
#PBS -l pmem=2g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)" 
conda activate rbase_env

cd /home/groups/nmrl/cov_analysis/downstream/
python update_heatmap_data.py ../mutation_files summary_file*