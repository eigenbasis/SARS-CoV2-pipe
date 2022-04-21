#!/bin/bash
#THE SCRIPT IS USED TO CHANGE STANDARD HEADERS IN CONSENSUSFIXER-PRODUCED FASTA FILE TO GISAID-ACCEPTED FASTA HEADERS.
consensus_path=${1} #PATH TO CONSENSUS.FASTA
sample_id_raw=${2} #SAMPLE ID
sample_id=$(python -c "print('${sample_id_raw}'.split('_',1)[1])")
#REGEX=$(CAT ${2}) #READING REGEX PATTERN FROM FILE
#SAMPLE_ID=$(ECHO $CONSENSUS_PATH | GREP -OP ${REGEX}) #EXTRACTING SAMPLE ID FROM PATH TO CONSENSUS.FASTA USING REGEX
sed -i "s/>Consensus_${sample_id_raw}_consensus_threshold_0.5_quality_30/>hCoV-19\/Latvia\/${sample_id}\/$(date +'%Y')/" ${consensus_path} #REPLACING CONSENSUS HEADER FOR GISAID FORMAT HEADER USING SED
