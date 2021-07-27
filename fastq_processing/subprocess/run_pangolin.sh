#!/bin/bash
ctrl_arg1=${1:-1} #if no parameter value provided, set 1st parameter value to 1 (perform Pangolin typing)
ctrl_arg2=${2:-0} #if no parameter value provided, set 2nd parameter value to 0 (skip Pangolin update)
eval "$(conda shell.bash hook)" #solves the issue of not being able to use conda activate <env-name> from bash script; ref: https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
conda activate pangolin # activating pangolin environment
now=$(date +"%m_%d_%Y")
cd ../reports/report_${now}
awk '{print}' *.fasta > ${now}_combined.fasta


if [ $ctrl_arg1 -eq 1 ] #if the first passed parameter is equal to 1, perform Pangolin typing
then
    if [ $ctrl_arg2 -eq 1 ] #if the 2nd passed parameter is equal to 1, perform Pangolin update
    then
        pangolin --update # checking for latest updates (slows the script down)
    elif [ $ctrl_arg2 -eq 0 ] #if the 2nd passed parameter is equal to 0, skip Pangolin update & report to the terminal that Pangolin update is skipped
    then 
    echo "Pangolin update option skipped."
    fi
    pangolin ${now}_combined.fasta --outfile ${now}_lineage_report.csv # running pangolin and providing time-dependent output file naming
elif [ $ctrl_arg1 -eq 0 ] #if the first passed parameter is equal to 0, skip Pangolin typing & report to the terminal that Pangolin typing is skipped
then
    echo "Pangolin typing skipped."
    if [ $ctrl_arg2 -eq 1 ]    
    then
        pangolin --update # checking for latest updates (slows the script down, uncomment once a day to check for updates)
    elif [ $ctrl_arg2 -eq 0 ] #if the 2nd passed parameter is equal to 0, skip Pangolin update & report to the terminal that Pangolin update is skipped
    then 
    echo "Pangolin update option skipped."
    fi
fi