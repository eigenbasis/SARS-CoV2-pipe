# SARS_CoV19_INHOUSE
Pipeline is used to generate .vcf files & human-readable csv files from Illumina PE150 fastq files of SARS-CoV19 WGS data.

## Tools & References
0. conda - https://docs.anaconda.com/
1. cutadapt 2.31 - https://doi.org/10.14806/ej.17.1.200
2. fastp 0.20.1 - https://doi.org/10.1093/bioinformatics/bty560
3. bwa 0.7.17-r1198-dirty - https://arxiv.org/abs/1303.3997
4. ivar 1.3.1 - https://doi.org/10.1186/s13059-018-1618-7
5. bedtools v.2.30.00 - https://doi.org/10.1093/bioinformatics/btq033
6. samtools 1.12 - https://doi.org/10.1093/bioinformatics/btp352
7. freebayes v0.9.21 - https://arxiv.org/abs/1207.3907
8. snpEff 5.0e - https://pcingola.github.io/SnpEff/adds/SnpEff_paper.pdf
9. Consensusfixer 0.4 - https://github.com/cbg-ethz/consensusfixer
10. abra 0.97 - https://doi.org/10.1093/bioinformatics/btu376
11. vcflib 1.0.2 - https://doi.org/10.1101/2021.05.21.445151 
12. Pangolin - https://doi.org/10.1038/s41564-020-0770-5
## Installation
### 1. Providing the folder structure
 - Clone the repository to the desired folder.
 - Create tools directory anywhere on your machine to store stand-alone tools.
 - In the **process_local_fastq.py** assign the path to the tools folder to the **tool_path** variable as string (e.g. tool_path = '/home/user/tools').
### 2. create cutadaptenv using conda
    conda create --name cutadaptenv
### 3. install dependencies under cutadaptenv
    conda activate cutadaptenv
#### python packages
	pip install pandas
	pip install scikit-allel
	pip install bokeh
#### cutadapt
	conda install -c bioconda cutadapt
#### fastp
	conda install -c bioconda fastp=0.20.1
#### bwa
 - [Installation](https://github.com/lh3/bwa.git)
 - [Adding to path](https://www.biostars.org/p/404164/)
#### ivar
	conda install -c bioconda ivar=1.3.1
#### bedtools
	conda install -c bioconda bedtools=2.30.00
#### samtools
	conda install -c bioconda samtools=1.12
#### vcflib
	conda install -c bioconda vcflib=1.0.2
#### freebayes
	conda install -c bioconda freebayes=0.9.21
### Stand-alone tools
***All stand-alone tools should be placed under tools directory created during step 1.***
#### snpEff 5.0e
 - [Installation](http://pcingola.github.io/SnpEff/download/)
 - [Create database](http://pcingola.github.io/SnpEff/se_buildingdb/)
 - [Reference genome](https://www.ncbi.nlm.nih.gov/nuccore/MN908947)
#### Consensusfixer 0.4
 - [Download](https://github.com/cbg-ethz/ConsensusFixer/releases/tag/0.4) the jar file to the tools directory .
#### abra 0.97
 - Create abra2 directory under tools directory.
 - [Download](https://github.com/mozack/abra/releases/tag/v0.97) jar file to the abra2 directory.
#### pangolin
 - [Install](https://cov-lineages.org/resources/pangolin/installation.html)
## Running the analysis
### Set the summary_path in subprocess/sample_stat_assembly.py script and plot_data in process_local_fastq.py script.
 - summary_path = f'path_to_base_stat_file/sample_summary_stats.csv'
 - plot_data = f'path_to_base_stat_file/sample_summary_stats.csv'
 - the sample_summary_stats.csv file must contain the following columns:
  	- receiving_lab_sample_id
 	- testing_lab
 	- normalized_sample_type
 	- sampling_date
 	- lineage
### Get a full path to the folder that contains fastq files to be analysed.
 - Fastq files should contain pair-end reads.
 - There should be distinct files for forward and reverse reads e.g. read_1 and read_2 files.
 - Fastq file name should contain sample id in 10-digit format that matches regex [0-9]{10}.
### Run process_local_fastq.py e.g. from terminal:
	python process_local_fastq.py -f path_to_fastq_containing_folder -m {metadata_file_name} -v {True/False}
 - metadata file name should be provided in order to include freshly-processed data into visuialization (file should be placed under resources/metadata)
 - If no metadata file name is provided, the pipeline will only generate mutation_report.csv (base statistics file will not be updated)
 - if v option is set to True - bokeh visualization will be displayed based on plot_data summary stats
### View output files in the reports folder
 - Mutation_report.csv file contains summary statistics about each sample.
    - Lineage and mutation columns report fraction of mutations that matched the given lineage/mutation filter for the given sample.
        - Filters are specified in fastq_processing/resources/report_filters.txt.
 - Ann.csv file for each sample contains list of detected mutations.
 - Depth_plot.html file contains coverage plot showing number of reads mapped to each position of reference genome.
### Metadata file
 - Metadata file should be placed under resources/metadata directory.
 - Values **must be stored as text**.
 - It should be a csv file containing the following columns:
 	- receiving_lab_sample_id
 		 - Must be a 10-digit code in that matches [0-9]{10} regex (e.g. 0000000001)
 	- testing_lab
 		 - Any text/numbers are accepted to represent testing laboratory.
 	- normalized_sample_type
 		 - Any text/numbers are accepted to represent sample types.
 	- sampling_date
 		 - must be in YYYY-MM-DD format (e.g. 2021-02-20)
