import json
import sys
import os

import numpy as np
import pandas as pd
#################################################################################################
#### README first! ##############################################################################
#################################################################################################
# Assume each matched protein contributes equally to the PSM of the peptide
# aka, moler concentration of each matched protein is equal to each other
# "Start Processing" section is the same like un-curated .py script except for PSM -> curated_PSM
# output
# json contains all the infos
# txt has the same content with Prof.YE code; \t seperated

print(f'Running {sys.argv[0]} ...')
print(f"Usage:\n\t$ {sys.argv[0]} ms_raw.txt")
#################################################################################################
#### Input ######################################################################################
#################################################################################################
# Get items from 559292
ID_S288c = '559292'
print(f'Organism taxid -> {ID_S288c}')

# IO
try:
    raw_data = sys.argv[1]
except:
    raw_data = '20250106_RVB2FLAG_raw.txt'    # for test
if raw_data == '':
    raise ValueError('No input files !')

try:
    f = open(raw_data)
except:
    print(f"Failed to open > {raw_data} <")

#################################################################################################
#### Curating PSMs ##############################################################################
#################################################################################################
# We do not modify the raw data file;

# Curate PSM of peptides mapping to more than 1 proteins
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

curated_PSM = dict()    # id to curated PSM
keys = peptide_chunk.keys()
for k in keys:
    matched_peptides = peptide_chunk[k][2:] # 1st row -> main line; 2nd row -> col head
    PSM_curated = 0
    for mp in matched_peptides:
        j = mp.split('\t')
        psm_peptide = int(j[3])     # PSM for peptide
        num_isoforms = int(j[5])    # number of matched proteins; for ids, see next column
        # Assume each matched protein contributes equally to the PSM of the peptide
        # aka, moler concentration of each matched protein is equal to each other
        PSM_curated += psm_peptide/num_isoforms

    curated_PSM[k] = PSM_curated

#################################################################################################
#### Start Processing ###########################################################################
#################################################################################################
# Read yeast_protein_table
table_path = os.path.join('..', '0_yeast_protein_table.json')
yeast_table = dict()
with open(table_path,'r', encoding='utf-8') as f:
    yeast_table = json.loads(f.read())

# Get main infos from raw MS data; ignore peptide items and modification infos
main_lines = []
with open(raw_data, 'r') as f:  # auto saved format is GBK
    main_lines = [l for l in f.readlines() if f"OX={ID_S288c}" in l]

# Generate MS table
ms_table = dict()
for i in main_lines:
    k = i.split('\t')
    uni_ID = k[0]
    # unique_peptides = k[5]
    peptides = int(k[6])
    #################################################################################################
    PSM = curated_PSM[uni_ID]     # Curation Here
    #################################################################################################
    MW = float(k[9])
    pI = float(k[10])

    ms_table[uni_ID] = yeast_table[uni_ID]
    ms_table[uni_ID]['Peptides'] = peptides
    ms_table[uni_ID]['PSM'] = PSM
    ms_table[uni_ID]['MW'] = MW
    ms_table[uni_ID]['pI'] = pI
    ms_table[uni_ID]['SCPHR'] = PSM/ms_table[uni_ID]['Residue_Number']*100

#################################################################################################
#### Cluster ranking by Cluster Average SCPHR ###################################################
#################################################################################################
ms_data = pd.DataFrame(ms_table).T
grouped = ms_data.groupby(by="Cluster_Name") # can be optimized in future

# Treat cluster as a whole: calculating cluster SCPHR
# sum(PSM of all subunits)/sum(aa length of all subunits)
cluster_PSM_sum = grouped["PSM"].sum()
cluster_PSM_sum.reset_index()
cluster_aa_sum = grouped['Residue_Number'].sum()
cluster_aa_sum.reset_index()
cluster_SCPHR = cluster_PSM_sum/cluster_aa_sum*100
cluster_SCPHR.reset_index()

# Write cluster_table: object is no longer single protein, but a complex
cluster_table = pd.DataFrame()
cluster_table['Cluster_Name'] = cluster_SCPHR.index
cluster_table['Cluster_PSM'] = cluster_PSM_sum.values
cluster_table['Cluster_Residue_Number'] = cluster_aa_sum.values
cluster_table['Cluster_SCPHR'] = cluster_SCPHR.values

# Sorted by SCPHR
cluster_table.sort_values(by='Cluster_SCPHR',ascending=False, inplace=True)
index = np.array(np.arange(0,len(cluster_table.index)))
cluster_table['Index'] = index
# print(cluster_table)

# Anno cluster index info to every protein item
map = pd.DataFrame()
map['Index'] = cluster_table['Index']
map.index = cluster_table['Cluster_Name']
map = map.to_dict()

keys = ms_table.keys()
for k in keys:
    try:
        ms_table[k]['Cluster_Index'] = map['Index'][ms_table[k]['Cluster_Name']]
        # ms_table[k]['Cluster_ind2name'] = {map['Index'][ms_table[k]['Cluster_Name']]: ms_table[k]['Cluster_Name']}
    except:
        ms_table[k]['Cluster_Index'] = 999

#################################################################################################
#### Cluster ranking by Cluster Total SCPHR ###################################################
#################################################################################################
ms_data = pd.DataFrame(ms_table).T
grouped = ms_data.groupby(by="Cluster_Name") # can be optimized in future

# Treat cluster as a whole: calculating cluster SCPHR
# sum(PSM of all subunits)/sum(aa length of all subunits)

cluster_SCPHR_total = grouped['SCPHR'].sum()
cluster_SCPHR_total.reset_index()

# Write cluster_table: object is no longer single protein, but a complex
cluster_table = pd.DataFrame()
cluster_table['Cluster_Name'] = cluster_SCPHR_total.index
cluster_table['Cluster_SCPHR_total'] = cluster_SCPHR_total.values

# Sorted by SCPHR
cluster_table.sort_values(by='Cluster_SCPHR_total',ascending=False, inplace=True)
index = np.array(np.arange(0,len(cluster_table.index)))
cluster_table['Index'] = index


# Anno cluster index info to every protein item
map = pd.DataFrame()
map['Index'] = cluster_table['Index']
map.index = cluster_table['Cluster_Name']
map = map.to_dict()

keys = ms_table.keys()
for k in keys:
    try:
        ms_table[k]['Cluster_Index_total'] = map['Index'][ms_table[k]['Cluster_Name']]
        # ms_table[k]['Cluster_ind2name_total'] = {map['Index'][ms_table[k]['Cluster_Name']]: ms_table[k]['Cluster_Name']}
    except:
        ms_table[k]['Cluster_Index_total'] = 999


#################################################################################################
#### Writing to figures #########################################################################
#################################################################################################
# Write ms_table to json file
to_json = input("Write annotated data to json file?:")
if to_json in ['y', 'yes', 'Y', 'Yes', 'YES']:
    output_file_name = raw_data[:-7] + 'ano.json'
    with open(output_file_name,'w', encoding='utf-8') as f:
        json.dump(ms_table, f,indent=2)
        print('ms table writed...')

to_txt = input("Write annotated data to txt file?:")
content = ''
if to_txt in ['y', 'yes', 'Y', 'Yes', 'YES']:
    # head is the same from Prof. YE
    content += 'Protein\tCluster\tAccession\tAnnotation\tMatched queries\tMatched peptide\tResidue Number\tSCPHR\tSubCluster\tDetails\t\n'
    keys = ms_table.keys()
    for k in keys:
        a = ms_table[k]["Gene_Name"]
        try:
            b = ms_table[k]["Cluster_Name"]
        except:
            b = ''
        c = ms_table[k]["Description"]
        d = ms_table[k]["PSM"]
        e = ms_table[k]['Peptides']
        f = ms_table[k]["Residue_Number"]
        g = ms_table[k]["SCPHR"]
        try:
            h = ms_table[k]["SubCluster_Name"]
            i = ms_table[k]["Details"]
        except:
            h = ''
            i = ''

        content += f'{a}\t{b}\t{k}\t{c}\t{d}\t{e}\t{f}\t{g}\t{h}\t{i}\t\n'

    file_name = raw_data[:-7] + 'ano.txt'
    with open(file_name,'w', encoding='utf-8') as f:
        f.write(content)
        print('ms table writed...')




