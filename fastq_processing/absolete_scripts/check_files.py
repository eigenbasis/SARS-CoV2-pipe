import os, sys, re
file_location = sys.argv[1]
#THE SCRIPT CAN BE USED TO CHECK IF THERE IS A MATCHING NUMBER OF VCF, BAM & FASTA FILES FOR ALL SAMPLE IDS CAPTURED BY REGEX

#GET REGEX FROM FILE
with open('./current_regex.txt', 'r') as regex_file:
    regex = regex_file.readline()

#INIT LIST TO STORE FILE NAMES FOR EACH FORMAT
fasta_list, vcf_list, bam_list = {}, {}, {}
for file in os.listdir(file_location):
    if os.path.isfile(f'{file_location}{file}') and ".fasta" in file and "combined" not in file:
        fasta_list[re.search(regex,file).group(0)] = f'{file_location}{file}'
    elif os.path.isfile(f'{file_location}{file}') and ".vcf" in file:
        vcf_list[re.search(regex,file).group(0)] = f'{file_location}{file}'
    elif os.path.isfile(f'{file_location}{file}') and ".bam" in file:
        bam_list[re.search(regex,file).group(0)] = f'{file_location}{file}'

#IF THERE IS EQUAL NUMBER OF FILES FOR EACH FORMAT
if len(fasta_list)==len(vcf_list)==len(bam_list):
    print('Equal number of all file types.')
else:
    print(len(fasta_list))
    print(len(vcf_list))
    print(len(bam_list))
    sys.exit('Mismatch in number of files of different types.')

#TEST ID MATCH BETWEEN FILE GROUPS (PRESUMPTIVE CHECKING - IF ID IS MATCHED IN FASTA-VCF-BAM THEN IT WILL MATCH BAM-VCF-FASTA AND VCF-BAM-FASTA)
for line in fasta_list.keys():
    try:
        vcf_list[line]
        bam_list[line]
    except KeyError:
        sys.exit(f'Fasta file id does not match vcf & bam id - {line}')
