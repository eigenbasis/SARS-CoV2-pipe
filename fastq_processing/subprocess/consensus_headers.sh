#!/bin/bash
#The script is used to change standard headers in Consensusfixer-produced fasta file to GISAID-accepted fasta headers.

consensus_path=${1} #path to consensus.fasta
sample_id=$(echo $consensus_path | grep -oP '[0-9]{10}') #extracting sample id from path to consensus.fasta using regex
sed -i "s/>CONSENSUS/>hCoV-19\/Latvia\/${sample_id}\/2021/" ${consensus_path} #replacing CONSENSUS header for GISAID format header using sed
