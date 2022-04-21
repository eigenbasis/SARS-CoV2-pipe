import pandas as pd, os, sys

#GET DICT MATCHING COV IDS AND REAL IDS
summary_df = pd.read_csv(sys.argv[1])
cov_real_map = {real:cov for real,cov in zip(list(summary_df['receiving_lab_sample_id']),list(summary_df['processing_id'])) if cov != "Z_BMC"}
# real_cov_map = {cov:real for cov,real in zip(list(summary_df['processing_id']),list(summary_df['receiving_lab_sample_id'])) if cov != "Z_BMC"}

#GET LIST OFD PATHS FOR FILES TO BE FIXED
fix_file = './files_to_fix.csv'
if os.path.isfile(fix_file):
    file_list = list(pd.read_csv(fix_file)['file_paths'])
else:
    file_list = []

    for (root,dirs,files) in os.walk(f'./failed/', topdown=True): 
        for name in files: 
            if 'COV' not in name: file_list.append(os.path.join(root, name))
    #SAVE LIST OF FILE PATHS TO FILE TO AVOID REPROCESSING WHOLE TREE EACH TIME

    result_df = pd.DataFrame.from_dict({'file_paths':file_list})
    result_df.to_csv('./files_to_fix.csv', header=True, index=False)



#ITERATE OVER FILE LIST
id_list = []
for file in file_list:
    #EXTRACT REAL ID
    file_name = file.split('/')[-1]
    #SPECIAL CASE - FILES STARTING WITH _
    if file_name[0] == "_":
        file_name = file_name[1:]
    if ".ann.csv" in file_name or ".vcf" in file_name:
        id = file_name.split(".")[0]
        # if id == '11':
        #     print(file_name)
    else:
        id = file_name.split("_")[0]
    id_list.append(id)
# id_df = pd.DataFrame.from_dict({'id_list':id_list})
# id_df.to_csv('./ids_to_fix.csv', header=True, index=False)
    #GET CORRESPONDING COV ID
    # missing_ids = []
    try:
        cov_id = cov_real_map[id]
        new_name = f'{os.path.dirname(file)}/{cov_id}_{file_name}'
    except:
        print(file)
    # CONSTRUCT NEW ID STRING

    # RENAME
    os.rename(file, new_name)
    os.system(f'chmod 775 {new_name}')