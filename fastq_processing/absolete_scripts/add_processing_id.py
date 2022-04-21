import os, sys, pandas as pd
#The script can be used to add processing id to historical files (e.g. COV_NNNNNN)

#PATHS
dir_path = sys.argv[1]
pair_list_path = sys.argv[2]

#PARSING INPUT
id_map = pd.read_csv(pair_list_path)
id_dict = dict(zip(id_map["real"], id_map['processing']))

name_dict = {}

#GET FILE NAME FROM VCF, BAM AND FASTA FILES (ERROR IF FILE IS IDENTIFIED BUT NOT FOUND IN ID_MAP TABLE)
for file in os.listdir(dir_path):
    if "vcf" in dir_path:
        try:
            # print(file, id_dict[file.split(".")[0]])
            name_dict[file] = f'{id_dict[file.split(".")[0]]}_{file}'
        except KeyError:
            print("ERROR"+file)
    elif "bam" in dir_path:
        try:
            # print(file, id_dict[file[:len(file)-11]])
            name_dict[file] = f'{id_dict[file[:len(file)-11]]}_{file}'
        except KeyError:
            print("ERROR"+file)
    elif "fasta" in dir_path:
        try:
            # print(file, id_dict[file[:len(file)-16]])
            name_dict[file] = f'{id_dict[file[:len(file)-16]]}_{file}'
        except KeyError:
            print("ERROR"+file)

#RENAMING FILES TO COVNNNNNN_ID FORMAT
for key in name_dict.keys():
    try :
        # print(f'{dir_path}{key}', f'{dir_path}{name_dict[key]}')
        os.rename(f'{dir_path}/{key}', f'{dir_path}/{name_dict[key]}')
        print(f"Source path ({key}) renamed to destination ({(name_dict[key])}) path successfully.")
  
    # If Source is a file 
    # but destination is a directory
    except IsADirectoryError:
        print(f"Source is a file ({key}) but destination ({(name_dict[key])}) is a directory.")
    
    # If source is a directory
    # but destination is a file
    except NotADirectoryError:
        print(f"Source is a directory ({key}) but destination ({(name_dict[key])}) is a file.")
    
    # For permission related errors
    except PermissionError:
        print(f"{key} - Operation not permitted.")
    
    # For other errors
    except OSError as error:
        print(f"{key} - {error}")
