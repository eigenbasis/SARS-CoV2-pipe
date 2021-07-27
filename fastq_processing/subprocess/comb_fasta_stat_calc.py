import re, pandas as pd

def readGenome(combined_fasta_path): # reading combined fasta file into genome_dict dictionary, mapping header (containing unique sample id) to the genome sequence string
    '''The function is used to read multifasta into dictionary - map consensus sequence to sample id.'''

    genome_dict = {} # initialize a dict to store header:sequence pairs
    with open(combined_fasta_path, 'r') as fasta_file: # open fasta file for reading
        for line in fasta_file: # go through the lines of the fasta file
            if line[0] == '>': # if line is a header
                fasta_header = line # store header until next header is reached
                genome_dict[fasta_header] = '' # initialize an empty sequence mapped to header
            if not line[0] == '>': # if line is a part of the genome sequence
                genome_dict[fasta_header] += line.rstrip() # add the line to the string that is mapped to the last reached header
    return genome_dict # dictionary containing header:genome pairs from combined fasta file


def get_nt_counting_stats(combined_fasta_path, skip=False, headers_changed = True):
    '''The function is used co calculate GC content and N content based on consensus sequence.'''

    if skip: # if skipping option is chosen
        print("Statistics calculation option skipped.")
        return None # return None for testing purposes
    elif not skip: # if skipping option is not chosen
        sequence_stats_df = pd.DataFrame(columns=["receiving_lab_sample_id","genome_length","genome_N_percentage","genome_GC_content"]) # initialize dataframe to store the statistics
        genome_dict_items = readGenome(combined_fasta_path).items() # get header:sequence pairs from readGenome() output
        
        for header, sequence in genome_dict_items: # looping through each header-sequence pair in genome dict
            sequence_length = len(sequence) # calculating genome length
            gc_content = round(100 * (sequence.count('G') + sequence.count('C')) / len(sequence), 2) # calculating gc-content for each sequence
            n_content = round(100 * (sequence.count('N') / len(sequence)), 2) # calculating invalid base content for each sequence
            if re.search(r'[A-Z]{2}[0-9]{4}\.B([0-9]{2}|[0-9])', header): new_header = re.search(r'[A-Z]{2}[0-9]{4}\.B([0-9]{2}|[0-9])', header).group(0)
            elif headers_changed: new_header = header.split("/")[2] # extracting sample id from fasta header to be used later in the column mapping process
            else: new_header = header.split("_")[1].strip() # extracting sample id from fasta header to be used later in the column mapping process if fasta headers were not changed
            sequence_stats_df = sequence_stats_df.append({"receiving_lab_sample_id":new_header, "genome_length":sequence_length, "genome_N_percentage":n_content, "genome_GC_content":gc_content}, ignore_index=True)
            #adding result row to the dataframe
        return sequence_stats_df # return resulting dataframe