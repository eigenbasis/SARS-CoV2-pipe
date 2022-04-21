regex=$(cat /home/groups/nmrl/cov_analysis/fastq_processing/current_regex.txt)
#THE SCRIPT CAN BE USED TO CHECK IF REGEX CAPTURES ALL SAMPLE IDS FOR FASTQ FILES UNDER SPECIFIED FOLDER AND SUBFOLDERS

#COUNT FASTQ FILES IN FOLDER & SUBFOLDERS
file_folder_path=${1}
echo $file_folder_path
echo $regex
fastq_count=$(ls -R ${file_folder_path} | grep fastq.gz | wc -l)
echo "There are ${fastq_count} fastq files in ${file_folder_path} folder."

#COUNT NUMBER OF REGEX CAPTURES FOR THE FASTQ FILES
detected_fq=$(ls -R ${file_folder_path} | grep fastq.gz | grep -P ${regex} | wc -l)
echo "${detected_fq} fastq file names fully/partially matched regex."

#COUNT NUMBER OF MISSES FOR FASTQ FILES
not_detected_fq=$(ls -R ${file_folder_path} | grep fastq.gz | grep -vP ${regex} | wc -l)
echo "${not_detected_fq} fastq file names did not match regex."

#PRINT LIST OF REGEX MATCHES TO STD-OUT
echo
echo "List of files that matched regex:"
matched_files=($(ls -R ${file_folder_path} | grep fastq.gz | grep -P ${regex}))
match_groups=($(ls -R ${file_folder_path} | grep fastq.gz | grep -oP ${regex}))
for i in "${!matched_files[@]}"; do
    printf "%s is in %s\n" "${match_groups[i]}" "${matched_files[i]}"
done
echo

#PRINT LIST OF REGEX MISSES TO STD-OUT
echo "List of files that did not match regex:"
unmatched_files=$(ls -R ${file_folder_path} | grep fastq.gz | grep -vP ${regex})
echo "${unmatched_files[*]}"
echo
