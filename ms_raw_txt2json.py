import json

#################################################################################################
#### Input ######################################################################################
#################################################################################################
# Get items from 559292
ID_S288c = '559292'
raw_data = '20250106_RVB2FLAG_raw.txt'

#################################################################################################
#### Convert ####################################################################################
#################################################################################################
peptide_chunk = dict()
with open(raw_data, 'r') as f:  # auto saved format is GBK
    lines = f.readlines()
    in_chunk =False
    id = ''
    for l in lines:
        if f"OX={ID_S288c}" in l:
            id =  l.split('\t')[0]
            peptide_chunk[id] = []
            in_chunk = True

        if in_chunk and (id != ''):
           peptide_chunk[id].append(l)
        else:
            id = ''
            in_chunk = True

# Write chunk
output_file_name = raw_data[:-4] + '.json'
with open(output_file_name, 'w', encoding='utf-8') as f:
    json.dump(peptide_chunk, f, indent=2)
    print("Raw data: txt -> json ...\nDone!")