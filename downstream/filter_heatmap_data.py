from sre_constants import OP_IGNORE
import pandas as pd
import numpy as np
import sys

#READING DATA & NUMBER OF CONSECUTIVE DATES TO USE IN FILTERING
data = pd.read_csv("n_by_date.csv")
unique_dates = sorted(data['sampling_date'].unique())
if len(sys.argv) < 2:
    sys.exit('Path to n_by_date.csv is required!')
consec_filter = int(sys.argv[1]) - 1

#GENERATING FILTER RANGE LISTS
filter_list = []
for i in range(0,len(unique_dates),consec_filter): #"SLIDING WINDOW" APPROACH - ITERATE OVER SORTED (ASC) LIST OF UNIQUE DATES
    if i+consec_filter < len(unique_dates): #IF LIST CONTAINS DATES TO INCLUDE INTO FILTER (TO AVOID OUT-OF-BOUNDS)
        filter_list.append([unique_dates[i+l] for l in range(consec_filter+1)]) #ADD CURRENT DATE AND NEXT n DATES (n = number of consecutive dates to observe mutations)
filter_array = np.array(filter_list, np.datetime64) #TO NUMPY ARRAY TO SPEED-UP FILTERING
ordered_blocks = [data[(data['sampling_date'].astype(np.datetime64).isin(filter_array[i]))] for i in range(len(filter_array))] #SPLIT RAW DATA INTO BLOCKS ORDERED BY CONSECUTIVE DATES

output_df = pd.DataFrame() #TO STORE RESULTS
for block in ordered_blocks: 
    mutation_filter = block[block['Fill'] == 0]['Mutation'].unique() #BLOCKS ARE DEFINED BY CONSECUTIVE DATES - IF ANY MUTATION HAS FREQUENCY 0, IT HAS TO BE EXCLUDED FROM THE GIVEN BLOCK
    filtered_blocks = block[~block['Mutation'].isin(mutation_filter)] #EXCLUDING MUTATIONS THAT WERE NOT OBSERVED IN A GIVEN BLOCK FOR n CONSECUTIVE DATES*
    #IF MUTATION OCCURS FOR n CONSECUTIVE DAYS IN ANOTHER BLOCK, IT WILL BE INCLUDED INTO THE RESULTING DATAFRAME
    output_df = output_df.append(filtered_blocks) #STORING FILTERING RESULTS

output_df.drop_duplicates(keep='last').to_csv('filtered_mutations.csv', header=True, index=False) #DUPLICATES ARISE FROM OVERLAP, WHICH IS NEEDED TO CHECK ALL CONSECUTIVE DATE