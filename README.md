# SARS_CoV19_INHOUSE
Pipeline is used to generate .vcf files & human-readable csv files from Illumina PE150 fastq files of SARS-CoV19 WGS data.

## Workflow summary
1. The analysis begins with **quality control of raw fastq files**:
   - Adapter sequences are removed using *cutadapt*[1] (screening for [default illumina adapter sequences](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm))
   - Reads are trimmed based on base quality (Phred score >= 30) and length (read length >= 50bp) using *fastp*[2]
2. **Reads that passed QC are aligned to the reference** using *bwa*[3].
   - The output from *bwa*[3] is passed to *samtools*[4] to produce binary alignment map (bam) file.
3. **Primer sequences are removed from bam file** using *ivar*[5].
   - Primer sequences should be provided in Browser Extensible Data (bed) format - generated 1 time for a primer set, from fasta file using *bwa*[3] & *bedtools*[6]. 
4. To improve alignment quality, **local realignment is performed** on primer-free bam file using *abra*[7].
   - bed file with realignment targets is created from primer-free bam file using *bedtools*[6].
5. **Alignment QC metrics are extracted** from raw and primer-free bam files using *samtools*[4].
6. **Variant-calling is performed** on post-realignment bam file using *freebayes*[8].
7. **Raw variants are filtered** using *vcflib/vcffilter*[9] based on quality ([QUAL](https://samtools.github.io/hts-specs/VCFv4.1.pdf) > 30) and sequencing depth ([DP](https://samtools.github.io/hts-specs/VCFv4.1.pdf) > 15).
8. **Filtered variants are annotated** using *snpEff*[10], using [genbank reference] (https://www.ncbi.nlm.nih.gov/nuccore/MN908947).
9. **Consensus sequence is generated** (in fasta format) from post-realignment bam file using *Consensufixer*[11].
    - invalid base (-) is called if coverage is less that 15.
10. **Sample id is added to fasta header and invalid bases are replaced with N** using bash scripts.
11. **Annotated vcf file is converted to csv format and coverage depth plot is generated**
from sequencing depth data extracted with *samtools* using inhouse-developed python scripts.
12. **Temporary files are deleted** after each sample is processed.
13. **Lineage assignment is performed** based on consensus sequence using *Pangolin*[12].
14. Inhouse-developed python & bash scripts are used to control the flow of analysis for multiple samples, generate summary report and visualize the results.
    

## Tools & References
1. cutadapt 2.31 - https://doi.org/10.14806/ej.17.1.200
2. fastp 0.20.1 - https://doi.org/10.1093/bioinformatics/bty560
3. bwa 0.7.17-r1198-dirty - https://arxiv.org/abs/1303.3997
4. samtools 1.12 - https://doi.org/10.1093/bioinformatics/btp352
5. ivar 1.3.1 - https://doi.org/10.1186/s13059-018-1618-7
6. bedtools v.2.30.00 - https://doi.org/10.1093/bioinformatics/btq033
7. abra 0.97 - https://doi.org/10.1093/bioinformatics/btu376 
8. freebayes v0.9.21 - https://arxiv.org/abs/1207.3907
9. vcflib 1.0.2 - https://doi.org/10.1101/2021.05.21.445151 
10. snpEff 5.0e - https://pcingola.github.io/SnpEff/adds/SnpEff_paper.pdf
11. Consensusfixer 0.4 - https://github.com/cbg-ethz/consensusfixer
12. Pangolin - https://doi.org/10.1038/s41564-020-0770-5
## Installation
### 1. Providing the folder structure
 - Clone the repository to your machine
 - Create tools directory anywhere on your machine to store stand-alone tools
 - In the **process_local_fastq.py** assign the path to the tools folder to the **tool_path** variable as string 
 	- e.g. **tool_path** = '/home/user/tools'
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
***All stand-alone tools should be placed under tools directory created in step 1.***
#### snpEff 5.0e
 - [Installation](http://pcingola.github.io/SnpEff/download/)
 - [Creating database](http://pcingola.github.io/SnpEff/se_buildingdb/)
 - [Reference genome](https://www.ncbi.nlm.nih.gov/nuccore/MN908947)
#### Consensusfixer 0.4
 - [Download](https://github.com/cbg-ethz/ConsensusFixer/releases/tag/0.4) the jar file & save it under **tools** directory
#### abra 0.97
 - Create abra2 directory under tools directory
 - [Download](https://github.com/mozack/abra/releases/tag/v0.97) the jar file & save it under **tools/abra2** directory
#### pangolin
 - [Installation](https://cov-lineages.org/resources/pangolin/installation.html)
## Running the analysis
### Set the paths to the **sample_summary_stats.csv** file
 - subprocess/sample_stat_assembly.py - **summary_path** = 'path_to_base_stat_file/sample_summary_stats.csv'
 - process_local_fastq.py - **plot_data** = 'path_to_base_stat_file/sample_summary_stats.csv'
 - The **sample_summary_stats.csv** file must contain the following columns (see **NMRL_SARS_CoV2_Inhouse_pipeline/test_summary_stats.csv** for example):
  	- receiving_lab_sample_id
  		- Must be a 10-digit code in that matches [0-9]{10} regex (e.g. 0000000001)
 	- testing_lab
 		- Any text/numbers are accepted
 	- normalized_sample_type
 		- Any text/numbers are accepted
 	- sampling_date
 		- Must be in *YYYY-MM-DD* format (e.g. 2021-02-20)
 	- lineage
 		- Any text/numbers are accepted
 - If any additional columns are available to display on dashboard table, please add:
   - **column_name** to the **metadata_list** variable in the **fastq_processing/plotting/lablin_data_serve.py** script.
   - Appropriate *name to be displayed on dashboard to the corresponding position of* **lv_name_list** in the **fastq_processing/plotting/lablin_data_serve.py** script.
   #### example:
   - **seq_date** column to be added from **sample_summary_stats.csv** to the dashboard table
   - **metadata_list** = ['receiving_lab_sample_id', 'testing_lab', 'normalized_sample_type', 'sampling_date', **'seq_date'**]
   - **lv_name_list** = ['Parauga ID', 'Testēšanas laboratorija', 'Par. veids', 'Par. ņemšanas datums', **'Sekv. datums'**]
### Get a full path to the folder that contains fastq files to be analysed.
 - Fastq files should contain pair-end reads
 - There should be distinct files for forward and reverse reads e.g. read_1 and read_2 files
 - Fastq file name should contain sample id in 10-digit format that matches regex [0-9]{10} (e.g. 1000000000)
### Run process_local_fastq.py e.g. from terminal:
	python process_local_fastq.py -f path_to_fastq_containing_folder -m {metadata_file_name} -v {True/False}
 - **metadata_file_name** should be provided in order to include freshly-processed data into visuialization 
 	- Metadata file should be placed under **fastq_processing/resources/metadata**
 - If no **metadata_file_name** is provided, the pipeline will only generate **mutation_report.csv** 
 	- **sample_summary_stats.csv** file will not be updated!
 - If *-v* option is set to *True* - bokeh visualization will be displayed based on **sample_summary_stats.csv** file (if data on 2 or more samples is available)
### View output files in the reports folder
 - **fastq_processing/reports/report_mmdd/mutation_report.csv** file contains summary statistics about each sample
    - *lineage* and *mutation columns* report *fraction of mutations that matched the given lineage/mutation filter* for the given sample
   	- Example: 
   	  - Filter *B.1.617 D614G,L452R,P681R* specifies 3 mutations that must be detected in order to assing B.1.617 lineage to the sample
   	  - If only 2 mutations are detected - the value in B.1.617 column in the **mutation_report.csv** will be around 0.67 (2/3 of the mutations were found)
    - Filters are defined in **fastq_processing/resources/report_filters.txt** file
   	- Filters can be added in the following format *Filter_name(space)Prot_mut_1,Prot_mut_2,...Prot_mut_N*
 - **reports/report_dd_mm_yyyy/sample_id_ann.csv** file for each sample contains list of detected mutations
 - **reports/report_dd_mm_yyyy/sample_id_depth_plot.html** file contains coverage plot showing number of reads mapped to each position of the reference genome
### Metadata file
 - Metadata file should be placed under resources/metadata directory
 - Values *must be stored as text*
 - It should be a csv file containing the following columns:
 	- **receiving_lab_sample_id**
 		 - Must be a 10-digit code in that matches [0-9]{10} regex (e.g. 0000000001)
 	- **testing_lab**
 		 - Any text/numbers are accepted
 	- **normalized_sample_type**
 		 - Any text/numbers are accepted
 	- **sampling_date**
 		 - must be in *YYYY-MM-DD* format (e.g. 2021-02-20)
 - See **resources/metadata/test_metadata.csv** for example
