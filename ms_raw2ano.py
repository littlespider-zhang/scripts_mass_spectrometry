import json
import sys
import os

import numpy as np
import pandas as pd

import csv
from collections import defaultdict

#################################################################################################
#### README first! ##############################################################################
#################################################################################################
# Be careful to protein isoforms, since some peptides can derive from different isoforms.
# If these peptides' hits are recorded,
#   *** WE ASSUME THE CONTRIBUTIONS OF ISOFOMRS ARE THE SAME. ***
# When calculating SCPHR, we usd curated PSMs (PSM/matched_proteins_number)
# e.g.  Peptide ELISNASDALDKIR can be assigned to Hsc82 and Hsp82(P15108;P02829), if its PSM is 1,
#       its contribution to Hsc82 is 0.5 (1/2), to Hsp82 is 0.5 (1/2)

class MassSpectrometryData():

    def __init__(self, filepath):
        # Preparation block
        self.yeast_protein_table = os.path.join('yeast_protein_table','0_yeast_protein_table.json')
        self.yeast_table = dict()
        self.taxid = '559292'  # S288c
        print(f'Organism taxid -> {self.taxid}')

        # Data block
        self.filepath = filepath
        self.curated_PSM = dict()
        self.ms_table_values = {"KEGG": '',
                                "Complex_Portal": '',
                                "Protein": '',
                                "UniProtID": '',
                                "Description": '',
                                "Structure": '',
                                "Matched_queries": 0,
                                "Matched_peptides": 0,
                                "Residue_Number": 0,
                                "SCPHR": 0.0,
                                "NAF": 0.0,
                                "SGD_ID": '',
                                "OLN": '',
                                }
        self.ms_table = dict()

        # Ribosomal subunits' isoform: generate from complex portal
        self.ribo_isoforms = [['O13516', 'P05755'], ['P02407', 'P14127'], ['P07280', 'P07281'],
                              ['P0C0W1', 'Q3E7Y3'], ['P0CX30'], ['P0CX31', 'P0CX32'], ['P0CX33'],
                              ['P0CX36', 'P0CX35'], ['P0CX38', 'P0CX37'], ['P0CX39', 'P0CX40'],
                              ['P0CX47', 'P0CX48'], ['P0CX51', 'P0CX52'], ['P0CX55', 'P0CX56'],
                              ['P23248', 'P33442'], ['P26786', 'P48164'], ['P32905', 'P46654'],
                              ['P35997', 'P38711'], ['P39516', 'P06367'], ['P39938', 'P39939'],
                              ['P41058', 'P41057'], ['Q08745', 'P46784'], ['Q3E754', 'P0C0V8'],
                              ['Q3E792', 'P0C0T4'], ['Q3E7X9', 'P0C0X0'], ['O14455', 'P05745'],
                              ['P02400', 'P05319'], ['P04449', 'P24000'], ['P05318', 'P10622'],
                              ['P05737', 'Q12213'], ['P05738', 'P51401'], ['P05740', 'P46990'],
                              ['P05743', 'P53221'], ['P05744', 'P41056'], ['P05748', 'P54780'],
                              ['P0C0W9', 'Q3E757'], ['P0C2H6', 'P0C2H7'], ['P0C2H8', 'P0C2H9'],
                              ['P0CH08'], ['P0CX23', 'P0CX24'], ['P0CX25'], ['P0CX27'],
                              ['P0CX41', 'P0CX42'], ['P0CX43', 'P0CX44'], ['P0CX46', 'P0CX45'],
                              ['P0CX50', 'P0CX49'], ['P0CX53', 'P0CX54'], ['P0CX82', 'P0CX83'],
                              ['P0CX84', 'P0CX85'], ['P0CX86'], ['P10664', 'P49626'],
                              ['P17076', 'P29453'], ['P26784', 'P26785'], ['P36105', 'P38754'],
                              ['P49166', 'P51402'], ['P56628', 'P05749'], ['P87262', 'P40525'],
                              ['Q02326', 'P05739'], ['Q02753', 'Q12672'], ['Q12690', 'P40212']]

        # eIFs and eEFs
        self.translation_factors = {"P02994": "eEF1α/EF-Tu",
                                    "P29547": "eEF1Bγ1",
                                    "P36008": "eEF1Bγ2",
                                    "P32324": "eEF2/EF-G",
                                    "P16521": "eEF3A",
                                    "P53235": "eIF2A",
                                    "P20459": "eIF2α",
                                    "P14741": "eIF2Bα",
                                    "P32502": "eIF2Bβ",
                                    "P12754": "eIF2Bδ",
                                    "P32501": "eIF2Bε",
                                    "P09064": "eIF2β",
                                    "P32481": "eIF2γ",
                                    "P38249": "eIF3a",
                                    "P06103": "eIF3b",
                                    "P32497": "eIF3c",
                                    "P10081": "eIF4A",
                                    "P07260": "eIF4E",
                                    "P39935": "eIF4G1",
                                    "P23301": "eIF5A1",
                                    "P39730": "eIF5B",
                                    "P12385": "eRF1",
                                    }

    def load_yeast_table(self):

        # Be careful to table path
        print(f"Protein table -> {self.yeast_protein_table}")
        with open(self.yeast_protein_table, 'r', encoding='utf-8') as f:
            self.yeast_protein_table = json.loads(f.read())

    def curate_PSM(self):
        #################################################################################################
        #### Curating PSMs ##############################################################################
        #################################################################################################
        # We do not modify the raw data file;
        # Curate PSM of peptides mapping to more than 1 proteins
        peptide_chunk = dict()
        with open(self.filepath, 'r') as f:  # auto saved format is GBK
            lines = f.readlines()
            in_chunk = False
            id = ''
            for l in lines:
                if (f"OX={self.taxid}" in l) or (f"OS=Saccharomyces cerevisiae" in l):
                    id = l.split('\t')[0]
                    peptide_chunk[id] = []
                    in_chunk = True

                if in_chunk and (id != ''):
                    peptide_chunk[id].append(l)
                else:
                    id = ''
                    in_chunk = True

        keys = peptide_chunk.keys()
        for k in keys:
            matched_peptides = peptide_chunk[k][2:]  # 1st row -> main line; 2nd row -> col head
            PSM_curated = 0
            for mp in matched_peptides:
                j = mp.split('\t')
                try:
                    psm_peptide = int(j[3])  # PSM for peptide
                    num_isoforms = int(j[5])  # number of matched proteins; for ids, see next column
                    # Assume each matched protein contributes equally to the PSM of the peptide
                    # aka, moler concentration of each matched protein is equal to each other
                    PSM_curated += psm_peptide / num_isoforms
                except:
                    print(f"\nHi: wrong here {mp}")
                    raise IOError("Empty lines in raw file!")

            self.curated_PSM[k] = PSM_curated

    def generate_ms_table(self):
        # Get main infos from raw MS data; ignore peptide items and modification infos
        main_lines = []
        with open(self.filepath, 'r') as f:  # auto saved format is GBK
            lines = f.readlines()
            main_lines = [l for l in lines if f"OX={self.taxid}" in l]
            if len(main_lines) == 0:
                print(f"\n\tOX={self.taxid} not found..\tSearching OS=Saccharomyces cerevisiae",end='')
                main_lines = [l for l in lines if f"OS=Saccharomyces cerevisiae" in l]
                if len(main_lines) != 0:
                    print("\tFound.")
                else:
                    raise IOError("No yeast proteins found. From generate_ms_table().")

        # Generate MS table
        for i in main_lines:
            k = i.split('\t')
            uni_ID = k[0]
            unique_peptides = k[5]
            peptides = int(k[6])
            #################################################################################################
            PSM = self.curated_PSM[uni_ID]  # Curation Here
            #################################################################################################
            MW = float(k[9])
            pI = float(k[10])
            try:
                self.ms_table[uni_ID] = self.yeast_protein_table[uni_ID]
            except:
                if '-' in uni_ID:
                    print(f"Pay attention:")
                    print(f"\t{uni_ID}",end='')
                    loc = uni_ID.find('-')
                    uni_ID = uni_ID[:loc]
                    self.ms_table[uni_ID] = self.yeast_protein_table[uni_ID]
                    print(f"\t-> {uni_ID}")


            self.ms_table[uni_ID]['Matched_queries'] = unique_peptides
            self.ms_table[uni_ID]['Matched_peptides'] = peptides
            self.ms_table[uni_ID]['PSM'] = PSM
            self.ms_table[uni_ID]['MW'] = MW
            self.ms_table[uni_ID]['pI'] = pI
            self.ms_table[uni_ID]['SCPHR'] = PSM / self.ms_table[uni_ID]['Residue_Number'] * 100

    def curate_gene_name(self):
        for uni_ID in self.ms_table.keys():

            # Ribosomal subunit proteins
            if (' ribosomal subunit protein' in self.ms_table[uni_ID]['Description']) and (
                    'RACK1' not in self.ms_table[uni_ID]['Description']):
                if self.ms_table[uni_ID]['Gene_Name'][-1] in ["A", "B"]:
                    self.ms_table[uni_ID]['Gene_Name'] = self.ms_table[uni_ID]['Description'].split()[-1][:-1]
                    self.ms_table[uni_ID]['Description'] = self.ms_table[uni_ID]['Description'][:-1]
                else:
                    self.ms_table[uni_ID]['Gene_Name'] = self.ms_table[uni_ID]['Description'].split()[-1]

            # Histones
            elif 'Histone H' in self.ms_table[uni_ID]['Description']:
                histone_anno = self.ms_table[uni_ID]['Description'].split()
                if len(histone_anno) == 2:
                    self.ms_table[uni_ID]['Gene_Name'] = histone_anno[1]

            # RNA Polymerases
            elif 'DNA-directed RNA polymerase' in self.ms_table[uni_ID]['Description']:
                RNAP_anno = self.ms_table[uni_ID]['Description'].split()[-1]
                self.ms_table[uni_ID]['Gene_Name'] = RNAP_anno[0] + RNAP_anno.lower()[1:]

            # Translation factors: eIF, eEF, eRF
            elif uni_ID in self.translation_factors.keys():
                self.ms_table[uni_ID]['Gene_Name'] = self.translation_factors[uni_ID]

            else:
                # Change CBF5 to Cbf5
                self.ms_table[uni_ID]['Gene_Name'] = self.ms_table[uni_ID]['Gene_Name'][0] + self.ms_table[uni_ID][
                                                                                                 'Gene_Name'].lower()[
                                                                                             1:]

    def combine_ribosomal_subunits(self):
        """
        e.g. for eL15A and eL15B, we will treat them as eL15
                UniProtID: a,b (string)
                Unique Peptides: a,b (string)
                Peptides: a,b (string)
                PSM: sum
                Residue_Number: average
                MW: a,b (string)
                pI: a,b (string)
                SCPHR = PSM/Residune_Number
        """
        for i in self.ribo_isoforms:
            if len(i) > 1:
                try:
                    check1 = self.ms_table[i[0]]
                    try:
                        check2 = self.ms_table[i[1]]
                        # print(f"Found A & B:\t{i[0]}, {i[1]}")
                        isoform_0 = self.ms_table.pop(i[0])
                        isoform_1 = self.ms_table.pop(i[1])
                        self.ms_table[f"{i[0]},{i[1]}"] = dict()
                        self.ms_table[f"{i[0]},{i[1]}"]['Description'] = isoform_0['Description']
                        self.ms_table[f"{i[0]},{i[1]}"]['OLN'] = f"{isoform_0['OLN']},{isoform_1['OLN']}"
                        self.ms_table[f"{i[0]},{i[1]}"][
                            'Structure'] = f"{isoform_0['Structure']},{isoform_1['Structure']}"
                        self.ms_table[f"{i[0]},{i[1]}"]['SGD_ID'] = f"{isoform_0['SGD_ID']},{isoform_1['SGD_ID']}"
                        self.ms_table[f"{i[0]},{i[1]}"]['Gene_Name'] = isoform_0['Description'].split()[-1]
                        self.ms_table[f"{i[0]},{i[1]}"]['Residue_Number'] = (int(isoform_0['Residue_Number']) + int(
                            isoform_1['Residue_Number'])) / 2
                        self.ms_table[f"{i[0]},{i[1]}"]['GO'] = f"{isoform_0['GO']}"
                        self.ms_table[f"{i[0]},{i[1]}"]['Cluster_KEGG'] = f"{isoform_0['Cluster_KEGG']}"
                        self.ms_table[f"{i[0]},{i[1]}"]['Cluster_Custom'] = f"{isoform_0['Cluster_Custom']}"
                        self.ms_table[f"{i[0]},{i[1]}"][
                            'Cluster_Complex_Portal'] = f"{isoform_0['Cluster_Complex_Portal']}"

                        self.ms_table[f"{i[0]},{i[1]}"][
                            'Matched_queries'] = f"{isoform_0['Matched_queries']},{isoform_1['Matched_queries']}"
                        self.ms_table[f"{i[0]},{i[1]}"][
                            'Matched_peptides'] = f"{isoform_0['Matched_peptides']},{isoform_1['Matched_peptides']}"

                        self.ms_table[f"{i[0]},{i[1]}"]['PSM'] = int(isoform_0['PSM']) + int(isoform_1['PSM'])
                        self.ms_table[f"{i[0]},{i[1]}"]['SCPHR'] = self.ms_table[f"{i[0]},{i[1]}"]['PSM'] / \
                                                                   self.ms_table[f"{i[0]},{i[1]}"][
                                                                       'Residue_Number'] * 100
                    except:
                        pass
                except:
                    pass

    def normalize(self,normalize2=''):
        # No normalize here
        if normalize2 == '':
            print("\nNo normalize.")
            for i in self.ms_table.keys():
                self.ms_table[i]['NAF'] = self.ms_table[i]['SCPHR']

    def write_annotated_data(self, custome_path=''):
        to_txt = input("Write annotated data to txt file?:")
        content = ''
        if to_txt in ['y', 'yes', 'Y', 'Yes', 'YES']:
            # Complex is from Cluster_custom
            content += 'Protein\tUniProtID\tComplex\tDescription\tStructure\tMatched queries\tMatched peptide\tResidue Number\tSCPHR\tComplex_Portal\tKEGG\tSGD_ID\tOLN\tNAF\n'
            keys = self.ms_table.keys()
            for k in keys:
                try:
                    a = self.ms_table[k]['Cluster_KEGG']
                except:
                    a = ''
                try:
                    b = self.ms_table[k]['Cluster_Complex_Portal']
                except:
                    b = ''
                try:
                    c = self.ms_table[k]['Gene_Name']
                except:
                    c = ''

                d = k
                try:
                    e = self.ms_table[k]['Description']
                except:
                    e = ''
                try:
                    f = self.ms_table[k]['Structure']
                except:
                    f = ''
                try:
                    g = self.ms_table[k]['Matched_queries']
                except:
                    g = ''
                try:
                    h = self.ms_table[k]['Matched_peptides']
                except:
                    h = ''
                try:
                    l = self.ms_table[k]['Residue_Number']
                except:
                    l = ''
                try:
                    m = self.ms_table[k]['SCPHR']
                except:
                    m = ''
                try:
                    n = self.ms_table[k]['SGD_ID']
                except:
                    n = ''
                try:
                    o = self.ms_table[k]['OLN']
                except:
                    o = ''
                try:
                    z = self.ms_table[k]['Cluster_Custom']
                except:
                    z = ''
                try:
                    y = self.ms_table[k]['NAF']
                except:
                    y = ''

                content += f'{c}\t{d}\t{z}\t{e}\t{f}\t{g}\t{h}\t{l}\t{m}\t{b}\t{a}\t{n}\t{o}\t{y}\n'

            file_name = custome_path
            if file_name == '':
                file_name = self.filepath.replace('raw', 'ano')

            # Do not use encoding='utf-8' here
            with open(file_name, 'w') as f:
                f.write(content)
                print(f'Writed to ***> {file_name} <***.')

    def annotate(self, output_path=''):
        self.load_yeast_table()

        print(f"Processing -> {self.filepath}",end='')
        # Do not change the sequence of following command
        self.curate_PSM()
        self.generate_ms_table()

        # Do not change the sequence of following command
        self.curate_gene_name()
        self.combine_ribosomal_subunits()

        # Do normalize ?
        self.normalize(normalize2='')   # '' means do not normalize

        print(f"Finised.\n\tTotal: {len(self.ms_table.keys())} proteins.")
        # write_data
        self.write_annotated_data(custome_path=output_path)

def get_raw_data_files(directory='data'):
    """Get list of .txt files with full paths"""
    return [os.path.join(directory, f)
            for f in os.listdir(directory)
            if f.lower().endswith('.txt') and f.lower().startswith('raw_')]

def get_ano_data_files(directory='data'):
    """Get list of .txt files with full paths"""
    return [os.path.join(directory, f)
            for f in os.listdir(directory)
            if f.lower().endswith('.txt') and f.lower().startswith('ano_')]

def process_all_raw_data():
    a = get_raw_data_files('data')
    for i in a:
        filepath = os.path.join(i)
        ms_data = MassSpectrometryData(filepath)
        ms_data.annotate()

def merge_spectral_files(input_files, output_file):
    """
    合并多个质谱数据文件的基础版本

    参数:
        input_files: 输入文件路径列表
        output_file: 输出文件路径
    """
    # 数据结构初始化
    protein_data = {}
    sample_order = []
    all_samples = set()

    # 处理每个输入文件
    for file_path in input_files:
        # 获取样本ID (假设文件名即样本ID)
        sample_id = os.path.splitext(os.path.basename(file_path))[0]
        if sample_id not in all_samples:
            sample_order.append(sample_id)
            all_samples.add(sample_id)

        with open(file_path, 'r') as f:
            # 读取表头
            headers = f.readline().strip().split('\t')

            # 获取列索引
            try:
                protein_idx = headers.index('Protein')
                complex_idx = headers.index('Complex')
                resnum_idx = headers.index('Residue Number')
                desc_idx = headers.index('Description')
                naf_idx = headers.index('NAF')
            except ValueError as e:
                print(f"文件 {file_path} 缺少必要列: {str(e)}")
                continue

            # 处理数据行
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < max(protein_idx, complex_idx, resnum_idx, desc_idx, naf_idx) + 1:
                    continue

                protein = parts[protein_idx]
                complex_val = parts[complex_idx]
                resnum = parts[resnum_idx]
                desc = parts[desc_idx]

                # 处理NAF值
                try:
                    naf = float(parts[naf_idx])
                except (ValueError, IndexError):
                    naf = 0.0

                # 初始化蛋白质记录
                if protein not in protein_data:
                    protein_data[protein] = {
                        'Group': complex_val,
                        'ResNum': resnum,
                        'SD': desc,
                        'samples': {s: 0.0 for s in sample_order}
                    }
                else:
                    # 保留首次出现的元数据
                    if not protein_data[protein]['Group']:
                        protein_data[protein]['Group'] = complex_val
                    if not protein_data[protein]['ResNum']:
                        protein_data[protein]['ResNum'] = resnum
                    if not protein_data[protein]['SD']:
                        protein_data[protein]['SD'] = desc

                # 更新当前样本的值
                protein_data[protein]['samples'][sample_id] = naf

    # 写入输出文件
    with open(output_file, 'w') as f:
        # 写表头
        header = ['Sample', 'Group', 'ResNum', 'NonzeroPuri', 'Average', 'SD'] + sample_order
        f.write('\t'.join(header) + '\n')

        # 写数据行
        for protein in sorted(protein_data.keys()):
            data = protein_data[protein]

            # 计算统计值 (占位符)
            samples = data['samples']
            nonzero = sum(1 for v in samples.values() if v > 0)
            avg = sum(samples.values()) / len(samples) if samples else 0

            # 准备行数据
            row = [
                protein,
                data['Group'],
                data['ResNum'],
                str(nonzero),
                f"{avg:.4f}",
                data['SD']
            ]

            # 添加样本值
            row += [f"{samples[s]:.4f}" if s in samples else '0.0000' for s in sample_order]

            f.write('\t'.join(row) + '\n')

def manual_order(custom,files,exclude):
    """See usage in ms_ano2matrix"""
    custom_order_1 = []
    for c in custom:
        for f in files:
            if c in f:
                custom_order_1.append(f)

    custom_order_2 = []
    for c in custom_order_1:
        for e in exclude:
            if e not in c:
                custom_order_2.append(c)

    return custom_order_2


if __name__ == "__main__":
    process_all_raw_data()
    # input_files = get_ano_data_files(directory='data')
    # merge_spectral_files(input_files, 'test.txt')
