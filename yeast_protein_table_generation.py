import json
import re
import requests

class KEGG_DataBase():
    def __init__(self):
        self.organism = 'sce'
        self.u2k=dict()     # UniProt ID -> KEGG Gene ID
        self.k2m=dict()     # KEGG ID -> Module ID
        self.m2n=dict()     # Module ID -> Module Name
        self.u2m=dict()   # UniProt ID -> Module ID

    def up2kegg(self):

        # UniProt ID  KEGG Gene ID
        print(f"Converting UniProt ID to KEGG ID for organism {self.organism} ...", end='')
        url = f"https://rest.kegg.jp/conv/{self.organism}/uniprot"
        lines = requests.get(url).text.split('\n')
        # Generate the id map
        for l in lines:
            if l != '':
                items = l.split('\t')
                up = items[0][3:].strip()
                name = items[1].strip()

                self.u2k[up] = name

        print(f"\tSuccess!")

        test = list(self.u2k.keys())[0]
        print(f"\teg: {test} -> {self.u2k[test]}")

    def kegg2module(self):

        # KEGG ID  Module ID
        print(f"Converting KEGG ID to Module ID for organism {self.organism} ...", end='')
        url = f"https://rest.kegg.jp/link/module/{self.organism}"
        lines = requests.get(url).text.split('\n')

        # Generate the id map
        for l in lines:
            if l != '':
                items = l.split('\t')
                up = items[0].strip()
                name = items[1][7:].strip()

                # Multiple assignments occurs
                try:
                    check = self.k2m[up]
                    self.k2m[up] += f",{name}"
                except:
                    self.k2m[up] = name

        print(f"\tSuccess!")

        test = list(self.k2m.keys())[0]
        print(f"\teg: {test} -> {self.k2m[test]}")

    def module2name(self):
        # Module ID  Module Name
        print(f"Converting Module ID to Module Name for organism ...", end='')
        url = f"https://rest.kegg.jp/list/module"
        lines = requests.get(url).text.split('\n')

        # Generate the id map
        for l in lines:
            if l != '':
                items = l.split('\t')
                up = items[0].strip()
                name = items[1].strip()

                # Multiple assignments occurs
                try:
                    check = self.m2n[up]
                    self.m2n[up] += f",{name}"
                except:
                    self.m2n[up] = name

        print(f"\tSuccess!")

        test = list(self.m2n.keys())[0]
        print(f"\teg: {test} -> {self.m2n[test]}")

    def up2mname(self):
        # Mapping UniProt ID to Module ID
        for up in self.u2k.keys():
            self.u2m[up] = dict()
            try:
                md = self.k2m[self.u2k[up]]
                mdname = ','.join([self.m2n[id.strip()] for id in md.split(',')])
            except:
                print(f"Error: {up} not assigned to KEGG modules ...")
                md = ''
                mdname = ''

            self.u2m[up]["Module_ID"] = md
            self.u2m[up]["Module_Name"] = mdname

    def generate(self):
        self.up2kegg()
        self.kegg2module()
        self.module2name()
        self.up2mname()

        # Writing
        write_or_not = input('Update protein_table now ?\n:')
        if write_or_not in ['y', 'yes', 'Y', 'Yes', 'YES']:
            with open(f"db_{self.organism}_KEGG_module_table.json", 'w', encoding='utf-8') as f:
                json.dump(self.u2m, f, indent=2)

class DataBases():

    def __init__(self):
        # S288c taxid
        self.ID = '559292'

        # Downloaded from https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz
        self.kb_complete_fasta = 'db_uniprot_sprot.fasta'
        self.length = 4316004  # recalculate this once database is updated
        # Calculated by
        # with open('test.txt','w', encoding='utf-8') as f:
        #     f.write(f"{length}")

        # Downloaded from https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/yeast.txt
        # using curl -o db_yeast_protein_infos_SGD.txt *url above*
        self.SGD_info = 'db_yeast_protein_infos_SGD.txt'

        # Downloaded from https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping_selected.tab.gz
        self.GO_info = 'db_YEAST_559292_idmapping_selected.tab'
        # PDB info is also here!!


        # Cluster maps
        # Download from https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv
        self.complex_portal = 'db_complex_portal_559292.tsv'
        # Generated from class KEGG_DataBase()
        self.KEGG = 'db_sce_KEGG_module_table.json'
        # self.cluster_map_file = '0_customized_protein_clusters_2500118.json'

        # empty table
        self.table = dict()

        # table infos
        self.table_values = {"Gene_Name": '',
                             "Description": '',
                             "OLN": '',
                             "SGD_ID": '',
                             "Structure": '',
                             "Residue_Number": 0,
                             "GO": '',
                             "Cluster_KEGG": '',
                             "Cluster_Complex_Portal": '',
                            }

        print('*' * 50)
        print(f"Makeing tables for all proteins of TaxID:{self.table}.")
        print('*' * 50)

    def update_db(self):
        pass

    def protein_list(self):
        """
            Get all proteins with given ID from UniProtKB, *ALL* proteins are here.
            Return a python list
            """
        proteins = ''
        with open(self.kb_complete_fasta, encoding='utf-8') as f:

            # iterate every line in the file
            for i in range(self.length):
                line = f.readline()
                loc1 = line.find('OX=')
                loc2 = line.find('GN=')
                if self.ID in line[loc1:loc2]:
                    proteins += f"{line.strip()}\n"

        return proteins.split('\n')

    def write_uniid_description(self):
        print(f"Writing IDs, Descriptions ...")
        protein_list = self.protein_list()
        for p in protein_list:
            # p should be like item below
            # >sp|Q04431|ZWINT_YEAST Outer kinetochore KNL1 complex subunit KRE28 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=KRE28 PE=1 SV=1
            x = p.split('|')
            try:
                uni_ID = x[1]  # uniprot id or uniprot accession get!
                OS_loc = x[2].find('OS=')
                y = x[2][:OS_loc].split(' ')
                description = ' '.join(y[1:]).strip()  # protein description get!
                self.table[uni_ID] = dict()
                self.table[uni_ID]['Description'] = description

                # yeast_table[uni_ID]['Gene_Name'] = ''
                # yeast_table[uni_ID]['SGS_ID'] = ''
                # yeast_table[uni_ID]['Residue Number'] = ''
                # yeast_table[uni_ID]['Structure'] = ''
            except:
                if p:
                    print(f"Something wrong in >>{p}<< ...")


        print("Done!")

    def write_name_oln_sgdid_structure_resinum(self):
        print("Writing Gene Name, OLN, SGD_ID, Structure Solved infos ..")
        # Make place
        current_ids = self.table.keys()
        for uni_ID in current_ids:
            self.table[uni_ID]['Gene_Name'] = ''
            self.table[uni_ID]['OLN'] = ''
            self.table[uni_ID]['SGD_ID'] = ''
            self.table[uni_ID]['Structure'] = ''
            self.table[uni_ID]['Residue_Number'] = 0
        
        
        # Start write infos
        start_at = 58   # Base on the file, becareful
        end_at = -5 # Base on the file, becareful
        with open(self.SGD_info, encoding='utf-8') as f:
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
                structure = j[-2]  # structure get! .. get?  -> BE CAREFUL ! see note about j
                if '(' in structure:
                    oln = j[-7]
                    uni_ID = j[-6]
                    SGD_ID = j[-4]  # SGD_ID get!
                    resi_num = j[-3]  # residue number get!
                else:
                    structure = '(no)'
                    oln = j[-6]
                    uni_ID = j[-5]
                    SGD_ID = j[-3]  # SGD_ID get!
                    resi_num = j[-2]  # residue number get!

                # Check if uni_ID is in self.table already
                if uni_ID not in current_ids:
                    self.table[uni_ID] = dict()
                    self.table[uni_ID] = self.table_values

                # Multi OLN for A protein exist
                if self.table[uni_ID]['OLN'] != '':
                        self.table[uni_ID]['OLN'] += f",{oln}"
                else:
                    self.table[uni_ID]['OLN'] = oln

                self.table[uni_ID]['Gene_Name'] = gene_name
                self.table[uni_ID]['SGD_ID'] = SGD_ID
                self.table[uni_ID]['Structure'] = structure
                self.table[uni_ID]['Residue_Number'] = int(resi_num)


        # Manually curated for TIF1, TEF1, EFT1 & HHF1
        # They are the same proteins with ***2 but located in other place in the genome
        print(f"Manually renaming:P10081，P02994，P32324，P02309")
        self.table['P10081']['Gene_Name'] = 'TIF1' # 'TIF1 or TIF2'
        self.table['P02994']['Gene_Name'] = 'TEF1' # 'TEF1 or TEF2'
        self.table['P32324']['Gene_Name'] = 'EFT1' # 'EFT1 or EFT2'
        self.table['P02309']['Gene_Name'] = 'HHF1' # 'HHF1 or HHF2'

        print("Done!")

    def write_GO(self):
        print(f"Writing GO infos ...")
        # Make place
        current_ids = self.table.keys()
        for uni_ID in self.table.keys():
            self.table[uni_ID]['GO'] = ''

        # Start writing infos
        with open(self.GO_info, encoding='utf-8') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                x = lines[i].split('\t')
                uni_ID = x[0]
                # print(x[5]) PDB id
                # GO ID list
                # GO:0005737; GO:0005829; GO:0005634; GO:0004032; GO:0004090; GO:0052588; GO:0045149; GO:0019568; GO:0034599; GO:0042843
                GO_ID_list = x[6]  # seperated by ;
                if 'GO' not in GO_ID_list:
                    GO_ID_list = ''

                # Check if uni_ID is in self.table already
                if uni_ID not in current_ids:
                    self.table[uni_ID] = dict()
                    self.table[uni_ID] = self.table_values


                self.table[uni_ID]['GO'] = GO_ID_list

        print('Done!')

    def write_KEGG(self):
        print("Writing KEGG module infos ..")
        # Make place
        current_ids = self.table.keys()
        for uni_ID in self.table.keys():
            self.table[uni_ID]['Cluster_KEGG'] = ''

        # Start writing infos
        with open(self.KEGG, 'r', encoding='utf-8') as f:
            KEGG_table = json.load(f)
            for uni_ID in KEGG_table.keys():

                # Check if uni_ID is in self.table already
                if uni_ID not in current_ids:
                    self.table[uni_ID] = dict()
                    self.table[uni_ID] = self.table_values

                self.table[uni_ID]['Cluster_KEGG'] = KEGG_table[uni_ID]['Module_Name']



        print('Done.')

    def write_Complex_Portal(self):
        print("Writing Complex_Portal infos ..")
        # Make place
        current_ids = self.table.keys()
        for uni_ID in self.table.keys():
            self.table[uni_ID]['Cluster_Complex_Portal'] = ''

        # Start writing infos
        with open(self.complex_portal, encoding='utf-8') as f:
            lines = f.readlines()
            for i in range(1,len(lines)):
                x = lines[i].split('\t')
                uni_IDs = [z.split('(')[0] for z in x[4].split(')|')]
                uni_IDs = [id for id in uni_IDs if "URS" not in id]
                uni_IDs = [id for id in uni_IDs if "CPX" not in id]
                uni_IDs = [id for id in uni_IDs if "CHEBI" not in id]
                uni_IDs = [id for id in uni_IDs if "EBI" not in id]
                current = len(uni_IDs)
                for index in range(current):
                    if '-PRO_' in uni_IDs[index]:
                        uni_IDs[index] = uni_IDs[index].split('-PRO_')[0]
                    if '[' in uni_IDs[index]:
                        loc_1 = uni_IDs[index].find('[')
                        loc_2 = uni_IDs[index].find(']')
                        z = uni_IDs[index][loc_1+1:loc_2].split(',')
                        uni_IDs[index] = z[0]
                        uni_IDs += z[1:]

                complex_name = x[1]

                for uni_ID in uni_IDs:
                    # Check if uni_ID is in self.table already
                    if uni_ID not in current_ids:
                        self.table[uni_ID] = dict()
                        self.table[uni_ID] = self.table_values

                    # print(f"uni_ID: {uni_ID}\t{complex_name}")
                    if self.table[uni_ID]['Cluster_Complex_Portal'] != '':
                        self.table[uni_ID]['Cluster_Complex_Portal'] += f",{complex_name}"
                    else:
                        self.table[uni_ID]['Cluster_Complex_Portal'] = complex_name

        print('Done!')

    def write_all(self):
        self.write_uniid_description()
        self.write_name_oln_sgdid_structure_resinum()
        self.write_GO()
        self.write_KEGG()
        self.write_Complex_Portal()

    def record_protein_isoforms(self):
        print("Recording protein isoforms from Complex_Portal infos ..")
        isoforms = ''
        with open(self.complex_portal, encoding='utf-8') as f:
            lines = f.readlines()
            for i in range(1, len(lines)):
                x = lines[i].split('\t')
                complex_name = x[1]

                uni_IDs = [z.split('(')[0] for z in x[4].split(')|')]
                uni_IDs = [id for id in uni_IDs if "URS" not in id]
                uni_IDs = [id for id in uni_IDs if "CPX" not in id]
                uni_IDs = [id for id in uni_IDs if "CHEBI" not in id]
                uni_IDs = [id for id in uni_IDs if "EBI" not in id]
                current = len(uni_IDs)
                for index in range(current):
                    if '-PRO_' in uni_IDs[index]:
                        uni_IDs[index] = uni_IDs[index].split('-PRO_')[0]

                    if '[' in uni_IDs[index]:
                        isoforms += f"{uni_IDs[index]}\t{complex_name}\n"

        with open('1_isoforms_from_complex_portal.txt', 'w', encoding='utf-8') as f:
            f.write(isoforms)

        print('Done!')

    def table2json(self,file_name="0_yeast_protein_table.json"):
        # Writing to json file
        write_or_not = input('Update protein_table now ?\n:')
        if write_or_not in ['y', 'yes', 'Y', 'Yes', 'YES']:
            with open(file_name, 'w', encoding='utf-8') as f:
                json.dump(self.table, f, indent=2)
                print(f"File writen in {file_name}")



if __name__ == '__main__':
    dbs = DataBases()
    dbs.write_all()
    dbs.table2json("0_yeast_protein_table.json")
