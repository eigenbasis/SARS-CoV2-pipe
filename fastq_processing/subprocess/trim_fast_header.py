import os, re, sys
from pathlib import Path
from datetime import date

full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[2]
dir_path = sys.argv[1]
count = 0
cur_year = date.today().year
for file in os.listdir(dir_path):
    if ".fasta" in file:
        with open(f'{dir_path}/{file}', 'r+') as fasta_file:
            fasta_content = fasta_file.read()
            name = file.split("/")[-1].split("_", 1)[1]
            if "_" in name:
                name = name[:len(name) - 16]
            print(name, file)
            if name is not None: count += 1
            # found = re.search('CONSENSUS', fasta_content)
            found = re.search(r'hCoV-19/Latvia/.*/[0-9]{4}',fasta_content)
            pattern = f'{name}'
            if found:
                renamed_content = re.sub(found.group(0), f'hCoV-19/Latvia/{name}/{cur_year}', fasta_content)
                print(renamed_content.split('\n')[0])
                print(count)
                with open(f'{dir_path}/{file}', 'w+') as output_file: output_file.write(renamed_content)
            elif not found:
                print(f'ID change failed: {file}')
                
