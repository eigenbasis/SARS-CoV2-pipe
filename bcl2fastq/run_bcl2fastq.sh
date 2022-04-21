#!/bin/sh
#PBS -N bcl2fastq
#PBS -l nodes=1:ppn=36,mem=216g
#PBS -l walltime=1:00:00
#PBS -q batch
#PBS -j oe
#PBS -A rakus

#This is a job script to be submitted though qsub.
#bcl2fastq command (line 22) can be edited according to Illumina documentation: 
#https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf
#N.B! Assign at least 32Gb of RAM for each core for optimal performance (-l nodes=1:ppn=12,pmem=3g or mem=36g or more if possible).

#Benchmarking:
#2021-12-08 - Oncoseq dataset (8.7Gb total size, 7 samples, about 2M Reads/sample - demultiplexing took about 1min 10 seconds with nodes=1:ppn=12,mem=36g)
#2021-12-09 - Covid dataset (51Gb total size, 196 samples, about 2M Reads/sample - demultiplexing took about 16min 30 seconds with nodes=1:ppn=24,mem=60g)
#2021-12-19 - Covid dataset (48Gb total size, 376 samples, about 1M Reads/sample - demultiplexing took about 40min with nodes=1:ppn=48:mem=256g; memory overhead - 128Gb was enough, cpu overhead - 24 was enough)

module load bcl2fastq/bcl2fastq-2.20.0

IN_DIR=$1
NAME=$2
SCRIPT_DIR=$3
mkdir -p $SCRIPT_DIR/demultiplexed/$NAME
FIN_DIR=$SCRIPT_DIR/demultiplexed/$NAME
date
bcl2fastq --barcode-mismatches 0 --no-lane-splitting -R $IN_DIR -o $FIN_DIR
date
chmod -R 775 $FIN_DIR
mkdir $FIN_DIR/Undetermined
mv $FIN_DIR/Undetermined_* $FIN_DIR/Undetermined/

###Rename Undetermined to avoid excessive processing of irrelevant data.
#for i in $(ls $OUT_DIR/Undetermined_S0_R*); do mv $OUT_DIR/Undetermined_S0_R1_001.fastq.gz $OUT_DIR/Undetermined_S0_R1_001.fastq1.gz; mv $OUT_DIR/Undetermined_S0_R2_001.fastq.gz $OUT_DIR/Undetermined_S0_R2_001.fastq1.gz ; done

###Move fastq files to the designated folder
#mv "$OUT_DIR/"*.fastq.gz $FIN_DIR/
