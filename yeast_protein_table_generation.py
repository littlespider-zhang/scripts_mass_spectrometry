import json

# Result is saved in 0_yeast_protein_table.json

# Extract protein list from uniprot_sprot database
# Downloaded from https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz
kb_complete_fasta = 'uniprot_sprot.fasta'
# Downloaded from https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/yeast.txt
# using curl -o yeast_protein_infos_SGD.txt *url above*
SGD_info = 'yeast_protein_infos_SGD.txt'
# Downloaded from https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping_selected.tab.gz
GO_info = 'YEAST_559292_idmapping_selected.tab'
# PDB info is also here!!
# S288c taxid
ID_S288c = '559292'
# Cluster map
cluster_map_file = '0_customized_protein_clusters_2500118.json'
# Create yeast table
yeast_table = dict()


print('*'*50)
print(f"Makeing tables for all proteins of TaxID:{ID_S288c}.")
print('*'*50)
#################################################################################################
#### Part 1: UniProt IDs & Descriptions #########################################################
#################################################################################################
length = 4316004    # recalculate this once database is updated
# Calculated by
# with open('test.txt','w', encoding='utf-8') as f:
#     f.write(f"{length}")
yeast_proteins = ''
with open(kb_complete_fasta, encoding='utf-8') as f:

    # iterate every line in the file
    for i in range(length):
        line = f.readline()
        loc1 = line.find('OX=')
        loc2 = line.find('GN=')
        if ID_S288c in line[loc1:loc2]:
            yeast_proteins += f"{line.strip()}\n"

print(f"Writing IDs, Descriptions ...")
protein_list = yeast_proteins.split('\n')
for p in protein_list:
    # p should be like item below
    # >sp|Q04431|ZWINT_YEAST Outer kinetochore KNL1 complex subunit KRE28 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=KRE28 PE=1 SV=1
    x = p.split('|')
    try:
        uni_ID = x[1]   # uniprot id or uniprot accession get!
        OS_loc = x[2].find('OS=')
        y = x[2][:OS_loc].split(' ')
        description = ' '.join(y[1:]).strip()   # protein description get!
        yeast_table[uni_ID] = dict()
        yeast_table[uni_ID]['Description'] = description

        # #
        # yeast_table[uni_ID]['Gene_Name'] = ''
        # yeast_table[uni_ID]['SGS_ID'] = ''
        # yeast_table[uni_ID]['Residue Number'] = ''
        # yeast_table[uni_ID]['Structure'] = ''
    except:
        if p:
            print(f"Something wrong in >>{p}<< ...")
print("Done!")
#################################################################################################
#### Part 2: Gene Name from SGD & SGD_ID & Residue Number & Structure Solved ####################
#################################################################################################
print("Writing Gene Name, SGD_ID, Residue Number, Structure Solved infos ..")
start_at = 58
end_at = -5
with open(SGD_info, encoding='utf-8') as f:
    infos = f.readlines()
    for i in infos[start_at:end_at]:
        # j should like be item below: (3) flag means 3D structure solved
        # ['AAC1', 'YMR056C', 'P04710', 'ADT1_YEAST', 'S000004660', '309', '13'] -> uni_id is at j[2]; without structure
        # OR ['AAC3', 'YBR085W', 'P18238', 'ADT3_YEAST', 'S000000289', '307', '(3)', '2'] -> with structure
        # OR ['ABF1;', 'BAF1;', 'OBF1;', 'REB2;', 'SBF1', 'YKL112W', 'P14164', 'ABF1_YEAST', 'S000001595', '731', '11'] -> id is at j[5]
        j = i.split()
        gene_name = j[0]  # gene name get! ... get ?
        if ';' in gene_name:
            gene_name = gene_name[:-1]

        # # Multiple names availble, BE CAREFUL; see note about j
        structure = j[-2]    # structure get! .. get?  -> BE CAREFUL ! see note about j
        if '(' in structure:
            uni_ID = j[-6]
            SGD_ID = j[-4]  # SGD_ID get!
            resi_num = j[-3]  # residue number get!
        else:
            structure = '(no)'
            uni_ID = j[-5]
            SGD_ID = j[-3]   # SGD_ID get!
            resi_num = j[-2] # residue number get!

        # For debugging
        #print([(gene_name, SGD_ID, uni_ID,resi_num, structure)])

        # Addind info to yeast table via UniProt ID
        try:
            yeast_table[uni_ID]['Gene_Name'] = gene_name
            yeast_table[uni_ID]['SGD_ID'] = SGD_ID
            yeast_table[uni_ID]['Residue_Number'] = int(resi_num)
            yeast_table[uni_ID]['Structure'] = structure
        except:
            print(f"\t{uni_ID} in {SGD_info} not found currently, added now...")
            yeast_table[uni_ID] = dict()
            yeast_table[uni_ID]['Note'] = f'No info in {SGD_info}\n'
print("Done!")
# Manually curated for TIF1, TEF1, EFT1 & HHF1
# They are the same proteins with ***2 but located in other place in the genome
yeast_table['P10081']['Gene_Name'] = 'TIF1 or TIF2'
yeast_table['P02994']['Gene_Name'] = 'TEF1 or TEF2'
yeast_table['P32324']['Gene_Name'] = 'EFT1 or EFT2'
yeast_table['P02309']['Gene_Name'] = 'HHF1 or HHF2'


#################################################################################################
#### Part 3: GO items ###########################################################################
#################################################################################################
print(f"Writing GO infos ...")
GO_info = 'YEAST_559292_idmapping_selected.tab'
with open(GO_info, encoding='utf-8') as f:
    lines = f.readlines()
    for i in range(len(lines)):
        x = lines[i].split('\t')
        uni_ID = x[0]
        # print(x[5]) PDB id
        # GO ID list
        # GO:0005737; GO:0005829; GO:0005634; GO:0004032; GO:0004090; GO:0052588; GO:0045149; GO:0019568; GO:0034599; GO:0042843
        GO_ID_list = x[6]   # seperated by ;
        if 'GO' not in GO_ID_list:
            GO_ID_list = ''

        try:
            yeast_table[uni_ID]['GO'] = GO_ID_list
        except:
            print(f"\t{uni_ID} in {GO_info} not found currently, added now...")
            yeast_table[uni_ID] = dict()
            yeast_table[uni_ID]['Note'] = f'No info in {GO_info}\n'
print('Done!')


#################################################################################################
#### Part 4: Cluster Infos ######################################################################
#################################################################################################
print("Writing customized cluster infos ...")
cluster_map = dict()
with open(cluster_map_file,'r', encoding='utf-8') as f:
    cluster_map = json.loads(f.read())

keys = yeast_table.keys()
for k in keys:
    try:
        gn = yeast_table[k]['Gene_Name']
    except:
        # print(f'No info in yeast_protein_table_for {k}')
        pass
    try:
        cluster_name = cluster_map[gn]['Cluster_Name']
        subcluster_name = cluster_map[gn]['SubCluster_Name']
        yeast_table[k]['Cluster_Name'] = cluster_name
        yeast_table[k]['SubCluster_Name'] = subcluster_name
    except:
        # print(f"No info in {cluster_map_file} for {k} ...")
        pass
    try:
        details = cluster_map[gn]['Details']
        yeast_table[k]['Details'] = details
    except:
        pass
print('Done!')



#################################################################################################
#### Finish! ####################################################################################
#################################################################################################
print('*'*50)
print(f"Yeast table has been made for taxid:{ID_S288c}.")
print(f"\t{len(yeast_table)} proteins in total !")
print('*'*50)


# Writing to json file
write_or_not = input('Update yeast_protein_table now ?\n:')
if write_or_not in ['y', 'yes', 'Y', 'Yes', 'YES']:
    with open("0_yeast_protein_table.json",'w', encoding='utf-8') as f:
        json.dump(yeast_table, f,indent=2)
        print(f"File writen in yeast_table.json")





