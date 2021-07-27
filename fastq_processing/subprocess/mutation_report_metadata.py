import sys, pandas as pd, re, time
ab_report_path = sys.argv[1]

def process_ab_metadata(data_path, output_path):
    '''The function is used to convert metadata AB report file into pipeline metadata csv file.'''

    timestr = time.strftime("%Y%m%d-%H%M%S")
    req_columns = {
        'Parauga numurs (Svītrkods)':'receiving_lab_sample_id',
        'Nosūtītājiestādes kods':'testing_lab',
        'Pacienta uzvārds':'testing_lab_sample_id',
        'Materiāla nosaukums':'sample_type',
        'Parauga noņemšanas datums':'sampling_date_1',
        'Parauga pieņemšanas laiks':'sampling_date_2',
        'Pacienta vecums':'age',
        'Pacienta dzimums V/S':'gender'
        }

    df = pd.read_excel(data_path)[list(req_columns.keys())].rename(columns=req_columns).fillna('0') #read data, keep only required columns, fill empty cells with 0
    df['testing_lab_sample_id'] = pd.to_numeric(df['testing_lab_sample_id'], errors='coerce').fillna('0') #removing non-numeric entries in testing_lab_sample_id field, fill empty cells with 0
    new_column = [df['sampling_date_1'][i] if df['sampling_date_1'][i] != '0' else df['sampling_date_2'][i] for i in range(len(df['receiving_lab_sample_id']))] #extracting sampling date from two columns provided

    for i, date_record in enumerate(new_column): 
        primary_split = str(date_record).split(" ")[0].strip(".") # for each value in the sampling date column - remove the time part
        if "." in primary_split: # if the date is formatted using dots (.)
            secondary_split = primary_split.split(".")[0:3] # get list of date components (date, month, year)
            if len(secondary_split[0]) == 2: # if the length of the first date component is 2
                new_column[i] = "-".join([secondary_split[2], secondary_split[1], secondary_split[0]]) # assume dd/mm/yyyy format - reformat to yyyy/mm/dd
            elif len(secondary_split[0]) == 4: # if the length of the first date component is 4
                new_column[i] = "-".join([secondary_split[0], secondary_split[1], secondary_split[2]]) # assume yyyy/mm/dd format - reformat to yyyy/mm/dd
        else: new_column[i] = primary_split

    for col in ['sampling_date_1', 'sampling_date_2']: del df[col] #replacing two sampling_date columns with currated sampling_date column
    df.insert(1, "sampling_date", new_column)

    regex_filters = {r"^(S|s)iekal":"Saliva", r".ztr":"Oro-pharyngeal swab",r"^(S|s)ekc":"Organ"} #regex filter dict to perform sample type normalization
    joint_regex = fr"({'|'.join(list(regex_filters.keys()))})" #used for primary match of sample type
    normalized_st_column = [] #initializing list to store the normalized sample type column

    for record in df["sample_type"]:
        if re.match(joint_regex, str(record)): # if record matches any or all regexes
            normalized = False
            for filter in regex_filters.keys(): #find the exact regex that it matches and normalize accordingly from regex_filter dict
                if re.match(filter, str(record)): 
                    normalized = True
                    normalized_st_column.append(regex_filters[filter])
                    break
            if not normalized: normalized_st_column.append("Other") #if sample type does not match any regex - normalize as type "Other"

    df.insert(4,"normalized_sample_type", normalized_st_column) # add normalized sample_type column to the all_excel_dump df
    df.testing_lab_sample_id = df.testing_lab_sample_id.astype(int).astype(str) #cast testing_lab_sample_id as integer and store as string
    df.to_csv(f'{output_path}/metadata_{timestr}.csv', index = False, header=True) #generate currated metadata file

    return f'{output_path}/metadata_{timestr}.csv'

process_ab_metadata(sys.argv[1], f'/mnt/c/cov_seq/scripts/fastq_processing/resources/metadata')